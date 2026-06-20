# step_04_build_protein_matrix.py の解説

## このコードの役割

sageの出力（ペプチドレベルの定量テーブル）をタンパク質レベルのマトリクスに集約し、後続の統計解析で使用可能な形式に変換するスクリプトです。FDRフィルタ、タンパク質推論、強度集約を行います。

## 全体フロー

```
ペプチドデータ(lfq.tsv) → FDRフィルタ → タンパク質推論 → 強度集約 → タンパク質マトリクス
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `results/sage_output/lfq.tsv` | sageによるペプチドレベル定量結果 |
| 入力 | `data/raw/human_proteome.fasta` | UniProt ID→Gene Symbol変換用 |
| 出力 | `results/protein_matrix_from_sage.csv` | タンパク質×サンプルマトリクス |

**所要時間の目安**: 約1-2分（ペプチド数による）

## セクション別解説

### セクション1: FASTAヘッダからGene Symbolマッピング構築

```python
def build_id_to_gene_map(fasta_path: str) -> dict:
    gn_re = re.compile(r"\bGN=([^\s]+)")
    # ヘッダ例: >sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens GN=TP53

    m = gn_re.search(desc_part)
    if m:
        gene = m.group(1)  # GN=TP53 から TP53 を抽出
    else:
        gene = entry_name.split("_")[0]  # P53_HUMAN から P53 を抽出
```

**目的**: UniProt accession（例：P04637）からGene Symbol（例：TP53）への変換マップを作成する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `re.compile(r"\bGN=([^\s]+)")` | GN=に続く遺伝子名を抽出する正規表現 |
| `header.split("|")` | UniProtヘッダの区切り文字でパースする |

**引数の意味**:
- `GN=`: Gene Name（遺伝子名）を示すUniProtヘッダの標準フィールド
- `\b`: 単語境界（GN=で始まる完全なフィールドをマッチ）
- `([^\s]+)`: 空白以外の文字列をキャプチャグループで取得

**つまずきやすいポイント**:
- 一部のタンパク質でGN=フィールドが欠落している場合がある
- Entry Name（P53_HUMAN）から遺伝子名（P53）を抽出する代替処理が必要

### セクション2: ペプチドのタンパク質推論

```python
def extract_primary_acc(proteins_field: str) -> str:
    # sageのproteinsカラム: "sp|ACC1|NAME;tr|ACC2|NAME;..."
    first = proteins_field.split(";")[0]  # 最初のタンパク質を選択
    parts = first.split("|")
    if len(parts) >= 2:
        return parts[1]  # accession部分を返す
```

**目的**: 複数のタンパク質にマッチするペプチドから代表タンパク質を選択する（Razor protein方式）。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `split(";")` | セミコロンで区切られた複数タンパク質を分割 |
| `split("|")` | UniProtフォーマット（sp|ACC|NAME）をパース |

**引数の意味**:
- Razor protein: 複数のタンパク質にマッチするペプチドを最も証拠の多いタンパク質に割り当てる手法
- 先頭選択: sageでは最もスコアの高いタンパク質が先頭に配置される

**つまずきやすいポイント**:
- 共有ペプチド（複数タンパク質に共通）の扱いが解析結果に影響する
- タンパク質推論の方針（Razor vs Protein groups）で結果が変わる

### セクション3: 強度集約とマトリクス構築

```python
# FDRフィルタ（q-value < 0.01）
df = df[df["q_value"] < Q_VALUE_THRESHOLD].copy()

# サンプル列を特定（.mzMLで終わる列）
sample_cols = [c for c in df.columns if c not in meta_cols]

# Gene Symbolごとにペプチド強度を合計
matrix = df.groupby("gene")[sample_cols].sum(min_count=1)

# 0値をNaNに変換（検出なし = 欠損値）
matrix = matrix.replace(0, pd.NA)
```

**目的**: 高品質なペプチドのみを選択し、Gene Symbolごとに強度を集約してタンパク質定量値を算出する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `df[df["q_value"] < 0.01]` | Q-value（FDR調整済みp値）でフィルタ |
| `groupby("gene").sum()` | 遺伝子ごとにペプチド強度を合計 |
| `min_count=1` | 最低1つの有効値がある場合のみ合計を計算 |
| `replace(0, pd.NA)` | 0値を明示的な欠損値に変換 |

**引数の意味**:
- `Q_VALUE_THRESHOLD = 0.01`: FDR < 1%の高信頼ペプチドのみ使用
- `sum()`: ペプチド強度を単純合算（MaxLFQ等の複雑な手法は使わない）
- `pd.NA`: pandasの標準的な欠損値表現

## 【深掘り】 タンパク質推論とは？

- **ひとことで**: ペプチド情報からタンパク質の存在と定量を推定する処理
- **定義**: 検出されたペプチドがどのタンパク質由来かを判定し、タンパク質レベルの定量値を算出
- **どんなとき使う**: 質量分析では直接タンパク質を測定できず、ペプチド経由でタンパク質を推定するため
- **落とし穴**: 共有ペプチド（複数タンパク質に共通）や類似配列による誤同定
- **参考**: 論文 Methods セクションの「Protein inference」

## 【深掘り】 Q-value (FDR) とは？

- **ひとことで**: 偽陽性率を制御するための統計的指標
- **定義**: False Discovery Rate（偽発見率）を各ペプチド同定に対して調整したp値
- **Q-value < 0.01**: 100個の同定のうち1個以下が偽陽性という意味
- **なぜ必要**: 大量のペプチド同定では偽陽性が蓄積するため、統計的制御が必須
- **参考**: 論文 Methods セクションの「Statistical analysis」

## よくある質問

**Q: なぜペプチドを合計してタンパク質定量にするの？**
A: 同じタンパク質由来のペプチドは発現量に比例するため、それらを合計することでタンパク質の総発現量を推定できます。より高度な手法もありますが、単純な合計でも有効です。

**Q: Q-value 0.01は厳しすぎませんか？**
A: プロテオミクスでは1%FDRが標準的です。より緩い閾値（5%等）も可能ですが、偽陽性が混入するリスクがあります。

**Q: Razor proteinとは何ですか？**
A: 複数のタンパク質にマッチするペプチドを、最も証拠の多い単一のタンパク質に割り当てる手法です。計算が単純で実用的なため広く使われています。

**Q: 自分のデータで試したい場合は？**
A: `LFQ_TSV`のパスを変更し、sage出力のlfq.tsvファイルを指定してください。FASTAファイルも対象種のプロテオームに変更が必要です。

## 出力されるマトリクスの形式

```
Protein,CRC001-N,CRC001-T,CRC002-N,CRC002-T,...
TP53,15234.2,45678.1,12345.6,39876.5,...
EGFR,89123.4,12456.7,76543.2,23456.8,...
MYC,5678.9,23456.1,8765.4,34567.2,...
```

- **行**: タンパク質（Gene Symbol）
- **列**: サンプル名（.mzML拡張子を除去）
- **値**: ペプチド強度の合計値
- **欠損**: 検出されなかったタンパク質はNaN

## 重要な注意点

- **データサイズ**: 数千タンパク質×数十サンプルのマトリクスが生成されます
- **メモリ使用量**: 大規模データセットでは数GB必要な場合があります
- **依存関係**: sage実行済みのlfq.tsvが前提
- **パフォーマンス**: ペプチド数が多い場合は処理時間が延びます

## 次のステップ

このスクリプトの出力（`protein_matrix_from_sage.csv`）は以下のステップで利用されます：
- Step 5: `step_04_preprocess.py`での前処理（ログ変換、欠損値補完）
- Step 6: `step_05_overview_visualization.py`での可視化
- Step 7: `step_06_differential_expression.py`での差分発現解析