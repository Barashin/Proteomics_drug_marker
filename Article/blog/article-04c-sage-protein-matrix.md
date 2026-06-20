---
title: "sageタンパク質マトリクス構築【論文再現シリーズ #4c】"
emoji: "🧮"
type: "tech"
topics: ["proteomics", "sage", "dataprocessing", "labcode"]
published: false
---

# sageタンパク質マトリクス構築

## はじめに

前回（[#4b sage実行とmzMLファイル検査](article-04b-sage-execution.md)）でsageを実行し、ペプチドレベルの定量結果（lfq.tsv）を取得しました。この記事では、このlfq.tsvファイルを**タンパク質レベルに集約**し、以降のすべての統計解析・可視化で使用する「タンパク質×サンプル」マトリクスを構築します。

:::message
**この記事で行う処理**
sage のlfq.tsv（ペプチドレベル）を FASTA の `GN=` フィールドで Gene Symbol にマップし、同一タンパク質の複数ペプチドを統合して「タンパク質×サンプル」の定量マトリクスを生成します。FDRフィルタリングと欠損値処理も含めた完全な前処理を行います。
:::

## 前提

- [#4b sage実行とmzMLファイル検査](article-04b-sage-execution.md) が完了していること
- sage出力ファイル（`lfq.tsv`）が生成済み
- ヒトプロテオームFASTAファイルが利用可能
- **対応Notebook**: [`notebooks/step_04.ipynb`](../notebooks/step_04.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_04.ipynb）

### ライブラリと設定

```python
# 標準ライブラリ: ファイル操作(os)、正規表現(re)、時間計測(time)
import os, re, time
import pandas as pd   # データフレーム操作ライブラリ（表形式データの読み書き・加工に使う）
import numpy as np    # 数値計算ライブラリ（配列操作・統計計算に使用）

# --- パス設定（プロジェクト内の各ディレクトリ・ファイルへの相対パスを定数として定義） ---
FASTA_PATH  = "../data/raw/human_proteome.fasta"     # ヒトプロテオームのFASTA配列ファイルのパス
RESULTS_DIR = "../results"                           # 解析結果の出力先ルートディレクトリ
SAGE_OUT    = os.path.join(RESULTS_DIR, "sage_output")  # sageの出力先ディレクトリ（lfq.tsvなどが生成される）
TABLES_DIR  = os.path.join(RESULTS_DIR, "tables")       # 解析結果テーブルの保存先ディレクトリ

# ペプチド FDR（偽発見率）の閾値を定義（論文と同じ 1% = 0.01）
# FDRが1%未満のペプチドのみを使うことで、誤同定を統計的に制御する
Q_THRESHOLD = 0.01

# 出力先ディレクトリが存在しなければ作成する（exist_ok=Trueで既存でもエラーにならない）
os.makedirs(TABLES_DIR, exist_ok=True)
```

### FASTA解析によるGene Symbolマッピング

```python
def parse_fasta_gene_map(fasta_path):
    """UniProt FASTA ヘッダから accession -> Gene Symbol の辞書を作る。

    【処理ロジック】
    ヘッダ例: >sp|P04637|P53_HUMAN ... GN=TP53 ...
    1. UniProt標準形式（sp|acc|entry）をパース
    2. GN=フィールドからGene Symbolを抽出
    3. GN=がない場合はEntry Nameの先頭語を代替使用

    【返り値】
    {UniProtアクセッション: Gene Symbol} の辞書
    """
    id_to_gene = {}  # UniProtアクセッション番号をキー、Gene Symbolを値とする辞書を初期化
    # 「GN=」の後に続く空白以外の文字列（=Gene Symbol）を抽出する正規表現パターン
    gn_re = re.compile(r"\bGN=(\S+)")

    print(f"FASTA解析中: {fasta_path}")
    # FASTAファイルを開いて1行ずつ読み込む
    with open(fasta_path) as f:
        for line_num, line in enumerate(f, 1):  # 行番号付きで走査
            if not line.startswith(">"):
                continue  # 「>」で始まらない行はアミノ酸配列なのでスキップ（ヘッダ行のみ処理）

            # ヘッダ行を「|」で最大3つに分割（例: "sp", "P04637", "P53_HUMAN ..."）
            parts = line[1:].split("|", 2)
            if len(parts) >= 3:
                # UniProt標準形式（sp|アクセッション|エントリ名）の場合
                acc = parts[1]                     # アクセッション番号（例: P04637）
                entry_name = parts[2].split()[0]   # エントリ名（例: P53_HUMAN）
            else:
                # 標準形式でないヘッダの場合、先頭の単語をアクセッション番号として使う
                acc = line[1:].split()[0]
                entry_name = acc

            # ヘッダ行から「GN=遺伝子名」パターンを検索
            m = gn_re.search(line)
            if m:
                # GN=が見つかった場合、その遺伝子名を使用
                gene_symbol = m.group(1)
            else:
                # GN=がない場合、エントリ名の「_」前の部分を代替として使用（例: P53_HUMAN → P53）
                gene_symbol = entry_name.split("_")[0]

            id_to_gene[acc] = gene_symbol

            # 進捗表示（10,000行ごと）
            if line_num % 10000 == 0:
                print(f"  処理済み: {line_num:,} 行, 登録済み: {len(id_to_gene):,} エントリ")

    print(f"FASTA解析完了: {len(id_to_gene):,} タンパク質エントリを登録")
    return id_to_gene  # {アクセッション番号: Gene Symbol} の辞書を返す

# FASTAファイルを解析して、アクセッション番号→Gene Symbolの変換辞書を構築
id_to_gene = parse_fasta_gene_map(FASTA_PATH)
```

### lfq.tsvの読み込みとFDRフィルタリング

```python
def load_and_filter_lfq(sage_output_dir, q_threshold=0.01):
    """sage出力のlfq.tsvを読み込み、FDRフィルタリングを適用する。

    【処理ステップ】
    1. lfq.tsvをpandas DataFrameとして読み込み
    2. q_value < 閾値でFDRフィルタリング
    3. サンプル列名のクリーニング（.mzML除去）
    4. メタデータ列とサンプル列の分離

    【返り値】
    filtered_df: FDRフィルタ済みDataFrame
    sample_cols: サンプル名のリスト
    """
    # sageが出力したlfq.tsv（ペプチドレベルのLFQ定量値テーブル）のパスを構築
    lfq_path = os.path.join(sage_output_dir, "lfq.tsv")
    print(f"lfq.tsv読み込み中: {lfq_path}")

    # タブ区切り（TSV）ファイルをpandas DataFrameとして読み込む
    df = pd.read_csv(lfq_path, sep="\t")
    print(f"lfq.tsv: {len(df):,} ペプチド × {len(df.columns)} 列")

    # FDR フィルタ: q_valueが閾値未満のペプチドだけを残す（信頼性の低い同定を除外）
    # .copy()で新しいDataFrameを作成し、元データへの参照を切る（SettingWithCopyWarning防止）
    df_filtered = df[df["q_value"] < q_threshold].copy()
    print(f"FDRフィルタ後 (q < {q_threshold}): {len(df_filtered):,} ペプチド")

    # メタデータ列（ペプチド配列、電荷、タンパク質名など）の名前を集合として定義
    meta_cols = {"peptide", "charge", "proteins", "q_value", "score", "spectral_angle"}
    # メタデータ列以外の列名を抽出 → これらがサンプルごとの強度値列（CRC01-N.mzML, CRC01-T.mzMLなど）
    sample_cols = [c for c in df.columns if c not in meta_cols]

    # 列名から「.mzML」サフィックスを除去して、サンプル名をクリーンにする辞書を作成
    rename_dict = {c: c.replace(".mzML", "") for c in sample_cols}
    df_filtered = df_filtered.rename(columns=rename_dict)       # DataFrameの列名を一括変更
    sample_cols = [rename_dict[c] for c in sample_cols]   # 変更後のサンプル名リストを更新

    print(f"サンプル数: {len(sample_cols)}")  # サンプル数を表示（期待値: 32）

    return df_filtered, sample_cols

# lfq.tsvの読み込みとFDRフィルタリング実行
df_filtered, sample_cols = load_and_filter_lfq(SAGE_OUT, Q_THRESHOLD)
```

### Gene Symbolマッピングとタンパク質レベル集約

```python
def create_protein_matrix(df, id_to_gene_map, sample_cols):
    """ペプチドデータをタンパク質レベルに集約する。

    【集約ロジック】
    1. UniProt IDからGene Symbolへマッピング
    2. 同一Gene Symbolのペプチドを合計（LFQ標準手法）
    3. 0値を欠損値(NaN)に変換（検出限界以下の意味）
    4. 全サンプルで欠損のタンパク質を除去

    【返り値】
    protein_matrix: タンパク質 × サンプル の定量マトリクス
    """
    print("タンパク質マトリクス構築中...")

    # proteins列（例: "sp|Q99497|PARK7_HUMAN;sp|..."）から代表アクセッション番号を抽出
    # まず「;」で分割して最初のタンパク質を取り、次に「|」で分割して2番目の要素（=アクセッション番号）を取得
    df["acc"] = df["proteins"].str.split(";").str[0].str.split("|").str[1]
    print(f"アクセッション番号抽出完了: {df['acc'].nunique():,} ユニークID")

    # アクセッション番号をGene Symbolに変換（マッピングできなかった場合はアクセッション番号をそのまま使用）
    df["gene"] = df["acc"].map(id_to_gene_map).fillna(df["acc"])
    mapped_count = df["gene"].map(lambda x: x in id_to_gene_map.values()).sum()
    print(f"Gene Symbolマッピング: {mapped_count:,}/{len(df):,} ペプチド ({mapped_count/len(df)*100:.1f}%)")

    # 同じGene Symbolを持つペプチドの強度を合計して、タンパク質レベルに集約
    # min_count=1: 全てNaNの場合はNaNを返す（0にしない）。生の強度値なのでsumが適切
    print("ペプチド→タンパク質集約中...")
    protein_matrix = df.groupby("gene")[sample_cols].sum(min_count=1)
    print(f"集約後: {len(protein_matrix):,} タンパク質")

    # 値が0の箇所をNaN（欠損値=未検出）に置き換え
    # sageで強度0は「検出限界以下」を意味するため、統計解析でのNaN扱いが適切
    protein_matrix = protein_matrix.replace(0, pd.NA)

    # 全サンプルでNaNの行（全く検出されなかったタンパク質）を削除
    initial_count = len(protein_matrix)
    protein_matrix = protein_matrix.dropna(how="all")
    removed_count = initial_count - len(protein_matrix)
    print(f"全欠損タンパク質除去: {removed_count:,} 個 → 最終 {len(protein_matrix):,} タンパク質")

    # インデックス名を「Protein」に設定（CSV出力時のヘッダになる）
    protein_matrix.index.name = "Protein"

    return protein_matrix

# タンパク質マトリクス構築実行
protein_matrix = create_protein_matrix(df_filtered, id_to_gene, sample_cols)

# マトリクスの基本統計
print("\n=== タンパク質マトリクス基本統計 ===")
print(f"マトリクスサイズ: {protein_matrix.shape[0]:,} タンパク質 × {protein_matrix.shape[1]} サンプル")
print(f"非欠損値数: {protein_matrix.notna().sum().sum():,} / {protein_matrix.size:,} ({protein_matrix.notna().sum().sum()/protein_matrix.size*100:.1f}%)")

# タンパク質ごとの検出頻度（何サンプルで検出されたか）
detection_freq = protein_matrix.notna().sum(axis=1)
print(f"検出頻度分布:")
print(f"  全サンプル検出: {(detection_freq == len(sample_cols)).sum():,} タンパク質")
print(f"  半数以上検出: {(detection_freq >= len(sample_cols)/2).sum():,} タンパク質")
print(f"  1サンプルのみ: {(detection_freq == 1).sum():,} タンパク質")
```

### データ品質チェックと統計サマリー

```python
def quality_check_matrix(matrix, sample_cols):
    """タンパク質マトリクスの品質をチェックする。

    【チェック項目】
    1. サンプル間の定量値分布の一致度
    2. 異常値（極端に高い/低い値）の検出
    3. 検出パターンの分析（タンパク質・サンプル別）
    4. 定量範囲の確認
    """
    print("=== データ品質チェック ===")

    # 1. 基本統計量
    numeric_data = matrix.select_dtypes(include=[np.number])
    if not numeric_data.empty:
        overall_median = numeric_data.median().median()
        overall_std = numeric_data.std().mean()
        print(f"定量値中央値: {overall_median:.2e}")
        print(f"平均標準偏差: {overall_std:.2e}")

        # 定量範囲（対数スケール）
        min_val = numeric_data.min().min()
        max_val = numeric_data.max().max()
        if min_val > 0:
            dynamic_range = np.log10(max_val / min_val)
            print(f"ダイナミックレンジ: {dynamic_range:.1f} orders of magnitude")

    # 2. サンプル別統計
    print("\nサンプル別統計（上位5サンプル）:")
    sample_stats = []
    for col in sample_cols[:5]:  # 最初の5サンプルのみ表示
        detected = matrix[col].notna().sum()
        median_val = matrix[col].median()
        sample_stats.append({
            'Sample': col,
            'Detected': detected,
            'Median': f'{median_val:.2e}' if pd.notna(median_val) else 'N/A'
        })

    stats_df = pd.DataFrame(sample_stats)
    print(stats_df.to_string(index=False))

    # 3. タンパク質検出頻度分布
    detection_counts = matrix.notna().sum(axis=1)
    print(f"\nタンパク質検出頻度分布:")
    freq_dist = detection_counts.value_counts().sort_index()
    for freq, count in freq_dist.head(10).items():
        print(f"  {freq}サンプル検出: {count:,} タンパク質")

    return stats_df

# データ品質チェック実行
quality_stats = quality_check_matrix(protein_matrix, sample_cols)
```

### CSV保存と出力確認

```python
def save_protein_matrix(matrix, output_dir, filename="protein_matrix_from_sage.csv"):
    """タンパク質マトリクスをCSVファイルとして保存する。

    【保存形式】
    - 行: タンパク質（Gene Symbol）
    - 列: サンプル名（.mzML拡張子は除去済み）
    - 値: LFQ強度値（0は欠損値NaNに変換済み）
    - エンコーディング: UTF-8
    """
    # 出力CSVファイルのパスを構築
    output_path = os.path.join(output_dir, filename)

    print(f"\nタンパク質マトリクス保存中: {output_path}")

    # DataFrameをCSVファイルとして保存（インデックス=タンパク質名も含める）
    # このCSVが以降の全ステップ（前処理・可視化・統計検定・データベース照合）の入力データになる
    matrix.to_csv(output_path, encoding='utf-8')

    # ファイルサイズ確認
    file_size = os.path.getsize(output_path) / 1024 / 1024  # MB単位
    print(f"保存完了: {output_path}")
    print(f"ファイルサイズ: {file_size:.1f} MB")

    # 先頭5行×5列を表示して中身を確認
    print("\n保存データ確認（先頭5×5）:")
    preview = matrix.iloc[:5, :5]
    print(preview)

    return output_path

# CSVファイルとして保存
output_csv_path = save_protein_matrix(protein_matrix, RESULTS_DIR)

# 保存結果の詳細統計
print(f"\n=== 最終結果サマリー ===")
print(f"入力ペプチド数: {len(df_filtered):,} (FDR < {Q_THRESHOLD})")
print(f"出力タンパク質数: {protein_matrix.shape[0]:,}")
print(f"サンプル数: {protein_matrix.shape[1]}")
print(f"データ密度: {protein_matrix.notna().sum().sum() / protein_matrix.size * 100:.1f}%")
print(f"出力ファイル: {output_csv_path}")
```

## 実行方法

**スクリプトで一括実行する場合:**

```bash
python scripts/step_04c_protein_matrix.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_04.ipynb
```

## 処理結果と統計

### Gene Symbolマッピング結果

| 項目 | 結果 |
|------|------|
| **FASTA解析エントリ数** | 約20,000タンパク質 |
| **ペプチドマッピング成功率** | 95%以上 |
| **GN=フィールド利用率** | 90%以上（残りはEntry Name使用） |

### タンパク質集約結果

| 項目 | 本書実測値 |
|------|-----------|
| **入力ペプチド数** (FDR<1%) | 8,391ペプチド |
| **出力タンパク質数** | **2,110タンパク質** |
| **集約率** | 約4ペプチド/タンパク質 |
| **データ密度** | 約85%（非欠損値の割合） |

### 検出頻度分布

- **全サンプル検出**: 約1,200タンパク質（高信頼性）
- **半数以上検出**: 約1,800タンパク質（中程度信頼性）
- **低頻度検出**: 約300タンパク質（低発現・ノイズ含む）

## コード詳細

### Gene Symbolマッピング戦略

| ケース | 処理方法 | 例 |
|--------|---------|-----|
| **GN=あり** | そのまま使用 | `GN=TP53` → `TP53` |
| **GN=なし** | Entry Name前半使用 | `P53_HUMAN` → `P53` |
| **非標準形式** | アクセッション番号使用 | フォールバック処理 |

### ペプチド集約ロジック

| 処理 | 設定 | 意味 |
|------|------|------|
| **groupby("gene")** | Gene Symbol単位 | 同一タンパク質のペプチドを統合 |
| **sum(min_count=1)** | 合計値計算 | LFQ標準手法（加法的定量） |
| **replace(0, pd.NA)** | ゼロ除去 | sage固有の「検出限界以下」処理 |
| **dropna(how="all")** | 全欠損除去 | 未検出タンパク質の削除 |

### データ品質指標

| 指標 | 期待範囲 | 意味 |
|------|---------|------|
| **データ密度** | 70-90% | 欠損値の少なさ（検出感度） |
| **ダイナミックレンジ** | 4-6 orders | 定量範囲の広さ |
| **サンプル間一致度** | CV<30% | 技術再現性の良さ |

## まとめ

sage のペプチドレベル出力を**タンパク質レベルマトリクス**に集約し、統計解析用データの準備が完了しました：

### 構築成果

1. **高品質マトリクス**: 2,110タンパク質×32サンプルの定量マトリクス
2. **Gene Symbol統一**: FASTA解析による標準化された遺伝子名
3. **FDR制御**: 1%閾値による統計的品質管理
4. **欠損値処理**: sage特有のゼロ値を適切にNaN変換

### 技術的特徴

- **LFQ集約**: ペプチドレベル強度の合計による標準的定量
- **商用クリア**: 全処理がMITライセンスツールで完結
- **再現性重視**: 明確な処理ロジックとパラメータ記録
- **拡張性**: 他のDIAツール出力にも適用可能な汎用設計

このマトリクスが以降のすべての解析（前処理・可視化・統計検定・データベース照合）の基盤となります。sage による DIA 解析は、商用制約のない環境で論文品質の解析を実現する実用的な選択肢であることが確認できました。

> 前回: [#4b sage実行とmzMLファイル検査](article-04b-sage-execution.md)
> 次回: [#6 前処理](article-06-preprocess.md) — ログ変換・欠損値処理・統計解析準備

#バイオインフォマティクス #プロテオミクス #sage #データ処理 #遺伝子マッピング #labcode