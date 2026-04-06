---
title: "商用利用OKのsage-proteomicsでDIA-MSデータを解析する【論文再現シリーズ #4】"
emoji: "🔮"
type: "tech"
topics: ["proteomics", "sage", "bioinformatics", "labcode"]
published: false
---

# 商用利用OKのsage-proteomicsでDIA-MSデータを解析する

## はじめに

この記事では、DIA-MSの生データからタンパク質を同定・定量する方法を解説します。論文（Toyota et al. 2025）では **DIA-NN v1.8.1** を使用していますが、DIA-NN は**商用利用に有料ライセンスが必要**です。

本シリーズでは **sage-proteomics（MIT ライセンス、完全無料・商用可）** で代替します。

| ツール | ライセンス | 商用利用 | 本書採用 |
|-------|-----------|---------|:---:|
| DIA-NN | 学術無料 / 商用有料 | ❌ | |
| FragPipe (MSFragger) | 学術無料 / 商用有料 | ❌ | |
| Spectronaut | 商用有料 | ⚠ | |
| **sage-proteomics** | **MIT** | **✅** | **✓** |
| OpenMS + AlphaPeptDeep | BSD + Apache 2.0 | ✅ | |
| MaxQuant / MaxDIA | 非商用無料 | ⚠ | |

## 前提

- [#1 環境構築](article-01-setup.md) が完了していること
- [#2 データ取得](article-02-data.md) で mzML ファイルが手元にあること

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
- **Rust 製**で DIA-NN 並みに高速（本書のデータ 18ファイル を **6.5分** で解析）
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

## sage の実行

```bash
sage \
  --fasta data/raw/human_proteome.fasta \
  --output_directory results/sage_output \
  --disable-telemetry-i-dont-want-to-improve-sage \
  scripts/sage_config.json \
  data/raw/raw_mzML/*.mzML
```

### 実行ログの要点

```
[INFO sage] generated 62769656 fragments, 3066221 peptides in 6029ms
[INFO sage] processing files 0 .. 5
[INFO sage] - search:     72700 ms (3507 spectra/s)
...
[INFO sage] finished in 262s
```

**約 6.5 分で 18 ファイルの DIA 解析が完了**します。DIA-NN と同等以上の速度です。

## 出力ファイル

```
results/sage_output/
├── results.json       (パラメータ記録)
├── results.sage.tsv   (PSM テーブル)
└── lfq.tsv            (ペプチド × サンプル定量マトリクス、27,586 行)
```

### lfq.tsv の中身

```
peptide   charge  proteins              q_value  score  CRC04-N.mzML   ...   CRC12-T.mzML
GLGTDEDTIIDIITHR  -1  sp|P08133|ANXA6_HUMAN  0.00122  0.815  94488098.3  ...   ...
DVVIC[+57.02146]PDASLEDAK  -1  sp|Q99497|PARK7_HUMAN  0.00207  0.742  28403038.6  ...   ...
```

- `charge = -1`: charge state は `combine_charge_states: true` により統合済み
- `q_value`: ペプチドレベル FDR。`< 0.01` でフィルタして論文と同じ 1% 基準に

## タンパク質マトリクスへの集約

sage の lfq.tsv はペプチドレベルなので、ダウンストリーム解析（第4章以降）で使うために **タンパク質レベルに集約** します：

```python
import re
import pandas as pd

df = pd.read_csv("results/sage_output/lfq.tsv", sep="\t")
df = df[df["q_value"] < 0.01]  # FDR 1% フィルタ

# UniProt ID から FASTA の GN= フィールドを使って Gene Symbol に変換
# （書籍側のコード参照: step_04_build_protein_matrix.py）

# Gene ごとに強度を合計
matrix = df.groupby("gene")[sample_cols].sum(min_count=1)
matrix.to_csv("results/protein_matrix_from_sage.csv")
```

## 結果サマリー（本書での実測値）

| 指標 | 本書 (sage, 18ファイル) | 論文 (DIA-NN, 32ファイル) |
|------|----------------------|------------------------|
| ペプチド数（生） | **27,586** | 記載なし |
| ペプチド数 (FDR<1%) | **8,391** | 記載なし |
| **タンパク質数** | **2,260** | **10,329** |
| 解析時間 | **6.5分** (Apple Silicon) | 記載なし |

タンパク質数は論文 (10,329) より少ない (2,260) ですが、これは以下の理由によります：

1. サンプル数 18 vs 32
2. sage は理論スペクトルベース（DIA-NN は深層学習ライブラリ予測）
3. sage は MBR 未実装

**商用完全クリアで論文と同等のパイプラインを学べる** のが本書の立ち位置です。

## まとめ

sage-proteomics は MIT ライセンスで、商用利用にも制約がない DIA 解析エンジンです。単一バイナリで導入が簡単、Rust 製で高速、library-free DIA をネイティブにサポートし、DIA-NN の商用代替として最適です。

> 前回: [#3 RAW→mzML変換](article-03-convert.md)
> 次回: [#5 前処理](article-05-preprocess.md) — Log2変換と欠損値補完

Sources:
- [sage GitHub](https://github.com/lazear/sage)
- [Lazear 2023 J. Proteome Res.](https://doi.org/10.1021/acs.jproteome.3c00486)
- [DIA-NN License](https://github.com/vdemichev/DiaNN)

#バイオインフォマティクス #プロテオミクス #sage #DIA #labcode
