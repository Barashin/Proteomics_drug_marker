---
title: "sage-proteomics基礎：DIA-MSの原理とセットアップ【論文再現シリーズ #4a】"
emoji: "🔮"
type: "tech"
topics: ["proteomics", "sage", "bioinformatics", "labcode"]
published: false
---

# sage-proteomics基礎：DIA-MSの原理とセットアップ

## はじめに

この記事では、DIA-MSの生データからタンパク質を同定・定量する方法を解説します。論文（Toyota et al. 2025）では **DIA-NN v1.8.1** を使用していますが、DIA-NN は**商用利用に有料ライセンスが必要**です。

本シリーズでは **sage-proteomics（MIT ライセンス、完全無料・商用可）** で代替します。

> **📝 INFO**
>
**この記事で行う処理**
DIA-MSの原理を理解し、sage-proteomicsのセットアップと設定を行います。DIA-MSデータをタンパク質定量マトリクスに変換するまでの全体像を学び、ヒトプロテオームのFASTAファイル準備からsage設定ファイル作成まで、実行準備を完了します。

## 生データから定量マトリクスまで：処理の全体像

この記事で行う処理を、「入力データ」「処理の中身」「出力データ」の3段階に分けて説明します。

### 入力：mzMLファイル（質量分析の生データ）とは

質量分析計に組織サンプルを入れると、サンプル中のペプチド（タンパク質の断片）が **質量/電荷比（m/z）** と **シグナル強度** のペアとして記録されます。これを「マススペクトル」と呼びます。

```
1つのmzMLファイルの中身（イメージ）:

スペクトル#1  → [(m/z=300.15, 強度=1200), (m/z=450.23, 強度=8500), ...]
スペクトル#2  → [(m/z=301.02, 強度=900),  (m/z=512.41, 強度=15000), ...]
...
スペクトル#3500 → [...]

→ 1ファイルあたり数千枚のスペクトル。
  「どのm/zにどれくらいのシグナルがあったか」の生記録。
```

この段階では **「どのタンパク質がどれくらいあるか」は全くわかりません**。m/zの数値の羅列があるだけです。DIA（Data Independent Acquisition）方式では、質量範囲を一定幅のウィンドウで網羅的にスキャンするため、サンプル中のほぼすべてのペプチドを漏れなく記録できるのが特徴です。

### 処理：sageが行う「同定」と「定量」

sage は以下の手順で、m/zの数値データを「タンパク質名 × 強度」に変換します。

**ステップ1：理論スペクトルの生成（in silico消化）**

ヒトの全タンパク質のアミノ酸配列が書かれたFASTAファイルを読み込み、トリプシン（タンパク質を切る酵素）で消化した場合に生じるペプチド断片を計算機上で生成します。さらに、各ペプチドが質量分析計の中で壊れたときに生じるフラグメントイオンの理論的なm/z値を計算します。

```
FASTAファイル（ヒト全タンパク質 約2万個の配列）
  ↓ トリプシン消化をシミュレーション
約300万個のペプチド断片（理論値）
  ↓ フラグメントイオンのm/zを計算
約6,300万個の理論フラグメント
```

**ステップ2：照合（同定）**

実測スペクトル（mzMLの中身）と理論スペクトルを比較し、「このm/zのピークはこのペプチドのフラグメントだ」とマッチングします。十分な数のフラグメントが一致すれば、そのペプチドが「同定」されます。

**ステップ3：定量（LFQ）**

同定されたペプチドについて、mzMLファイル内のシグナル強度（ピーク面積）を計測し、各サンプルにおける量を数値化します。これが **Label-Free Quantification（LFQ）** です。

### 出力：タンパク質×サンプルの定量マトリクス

sageの直接の出力はペプチドレベルの定量値（`lfq.tsv`）です。これを後段のスクリプトでタンパク質レベルに集約し、最終的に以下の形のCSVファイルが得られます。

```
               CRC01-N    CRC01-T    CRC02-N    CRC02-T   ...  CRC16-T
タンパク質A     18.5       22.3       19.1       21.8      ...  23.0
タンパク質B     15.2       14.8       15.9       15.1      ...  14.5
タンパク質C      NaN       12.3       11.8       13.0      ...  12.7
...
（2,110行 × 32列）
```

- **行**：同定されたタンパク質（2,110個）
- **列**：各サンプル（16患者 × Normal/Tumor = 32サンプル）
- **値**：そのサンプルにおけるタンパク質の強度（量）。NaNは検出限界以下で検出されなかったことを意味します

このマトリクスが以降すべてのステップ（前処理→可視化→統計検定→データベース照合）の入力データになります。

---

| ツール | ライセンス | 商用利用 | 本書採用 |
|-------|-----------|---------|:---:|
| DIA-NN | 学術無料 / 商用有料 | ❌ | |
| FragPipe (MSFragger) | 学術無料 / 商用有料 | ❌ | |
| Spectronaut | 商用有料 | ⚠ | |
| **sage-proteomics** | **MIT** | **✅** | **✓** |
| OpenMS + AlphaPeptDeep | BSD + Apache 2.0 | ✅ | 後述 |
| MaxQuant / MaxDIA | 非商用無料 | ⚠ | |

## 前提

- [#1 環境構築](article-01-setup.md) が完了していること
- [#2 データ取得](article-02-data.md) で mzML ファイルが手元にあること
- ヒトプロテオームのFASTAファイル（下記「FASTAファイルの準備」で取得）
- **対応Notebook**: [`notebooks/step_04.ipynb`](../notebooks/step_04.ipynb) — この記事のコードをセルごとに実行できます

## FASTAファイルの準備（ヒトプロテオーム配列）

sage はmzMLファイルだけでは動きません。「どのタンパク質が存在しうるか」の候補リストとして、**ヒト全タンパク質のアミノ酸配列を記録したFASTAファイル**が必要です。

### FASTAファイルとは

FASTAファイルは、タンパク質（または遺伝子）の配列をテキストで記述するフォーマットです。`>` で始まるヘッダ行と、それに続くアミノ酸配列の組で構成されます。

```
>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens GN=TP53
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLPGRNSFEV
...
```

- `sp|P04637|P53_HUMAN`: UniProt ID とエントリ名
- `GN=TP53`: 遺伝子シンボル（Gene Symbol）。sageの出力をタンパク質名に変換する際に使います
- 以降の行: アミノ酸の1文字表記による配列

ヒトの場合、約 **20,000 個** のタンパク質エントリが含まれます。

### UniProtからダウンロードする

FASTAファイルは **UniProt**（https://www.uniprot.org）からダウンロードできます。ヒトの参照プロテオーム（reviewed、Swiss-Prot）を取得します。

**ブラウザの場合:**

1. https://www.uniprot.org/uniprotkb?query=reviewed:true+AND+organism_id:9606 にアクセス
2. ページ上部の **「Download」** ボタンをクリック
3. Format で **「FASTA (canonical)」** を選択してダウンロード

**コマンドラインの場合:**

```bash
# UniProt からヒト reviewed プロテオームをダウンロード
curl -o data/raw/human_proteome.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=reviewed:true+AND+organism_id:9606"

# 確認: エントリ数を表示
grep -c "^>" data/raw/human_proteome.fasta
# → 約 20,000 エントリ
```

> **⚠️ WARNING**
>
**FASTAファイルのバージョンについて**: UniProtは定期的に更新されるため、ダウンロード時期によってエントリ数が若干異なります。論文と完全に同じ結果を再現したい場合は、論文で使用されたバージョンのFASTAを使うのが理想ですが、本シリーズでは最新版で問題ありません。


## sage-proteomics とは

| 項目 | 内容 |
|------|------|
| 開発者 | Michael Lazear |
| 初版 | 2022年 |
| 言語 | Rust（単一バイナリ、依存なし） |
| ライセンス | **MIT**（商用完全OK） |
| DIA サポート | ✅ library-free DIA ネイティブ対応 |
| 論文 | Lazear MR. *J. Proteome Res.* 2023. DOI: 10.1021/acs.jproteome.3c00486 |
| GitHub | https://github.com/lazear/sage |

**sage の強み**:

- **単一バイナリ**でインストールが超簡単
- **Rust 製**で DIA-NN 並みに高速（本書のデータ 32ファイル を **10.7分** で解析）
- **MIT ライセンス**なので企業・出版物でも自由に使える
- **library-free** で FASTA を渡すだけで in silico 消化から定量まで完結

## sage のインストール

```bash
micromamba activate crc-proteomics
micromamba install -c bioconda -c conda-forge sage-proteomics -y
sage --version
# → sage 0.14.6
```

## sage 設定ファイル (sage_config.json)

論文の DIA-NN パラメータに合わせた設定ファイルを作ります。

```json
{
  "database": {
    "bucket_size": 32768,
    "fragment_min_mz": 200.0,
    "fragment_max_mz": 1800.0,
    "peptide_min_mass": 600.0,
    "peptide_max_mass": 4000.0,
    "enzyme": {
      "missed_cleavages": 1,
      "min_len": 7, "max_len": 45,
      "cleave_at": "KR", "restrict": "P", "c_terminal": true
    },
    "static_mods": { "C": 57.02146 },
    "variable_mods": {},
    "max_variable_mods": 1,
    "ion_kinds": ["b", "y"],
    "decoy_tag": "rev_", "generate_decoys": true,
    "fasta": "data/raw/human_proteome.fasta"
  },
  "quant": {
    "lfq": true,
    "lfq_settings": {
      "peak_scoring": "Hybrid",
      "integration": "Sum",
      "spectral_angle": 0.7,
      "ppm_tolerance": 10.0,
      "combine_charge_states": true
    }
  },
  "precursor_tol": { "ppm": [-10, 10] },
  "fragment_tol":  { "ppm": [-10, 10] },
  "precursor_charge": [2, 4],
  "isotope_errors": [0, 0],
  "deisotope": false,
  "chimera": true,
  "wide_window": true,
  "predict_rt": false,
  "min_peaks": 15, "max_peaks": 150,
  "min_matched_peaks": 4,
  "max_fragment_charge": 2,
  "report_psms": 1,
  "output_directory": "results/sage_output"
}
```

### 論文の DIA-NN パラメータとの対応

| 論文（DIA-NN） | sage の対応 |
|---------------|------------|
| `--fasta-search`（library-free） | `wide_window: true` + FASTA指定 |
| `--cut K*,R*`（トリプシン） | `enzyme.cleave_at: "KR"` |
| `--missed-cleavages 1` | `enzyme.missed_cleavages: 1` |
| `--min-pep-len 7` | `enzyme.min_len: 7` |
| `--max-pep-len 45` | `enzyme.max_len: 45` |
| `--pr-charges 2-4` | `precursor_charge: [2, 4]` |
| `--min-fr-mz 200` / `--max-fr-mz 1800` | `fragment_min_mz: 200` / `fragment_max_mz: 1800` |
| MS1/MS2 accuracy 10 ppm | `precursor_tol.ppm: [-10,10]` / `fragment_tol.ppm` |
| Cys carbamidomethylation | `static_mods: { "C": 57.02146 }` |
| FDR <1% | 出力の `peptide_q < 0.01` で後段フィルタ |
| MBR (Match Between Runs) | sage 0.14 では実装なし（将来実装予定） |

## まとめ

DIA-MSプロテオミクス解析の原理を理解し、sage-proteomicsのセットアップが完了しました。FASTAファイルの準備から設定ファイルの作成まで、実行に必要な環境が整いました。

> 前回: [#3 データ変換](article-03-convert.md)
> 次回: [#4b sage実行と結果解析](article-04b-sage-execution.md) — sageの実行とタンパク質マトリクス生成

#バイオインフォマティクス #プロテオミクス #sage #DIA-MS #labcode