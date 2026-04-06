---
title: "無料ツールだけでDIA-MSプロテオミクス環境を構築する【論文再現シリーズ #1】"
emoji: "🛠"
type: "tech"
topics: ["proteomics", "python", "bioinformatics", "labcode"]
published: false
---

# 無料ツールだけでDIA-MSプロテオミクス環境を構築する

## はじめに

この記事では、DIA-MSプロテオミクス解析に必要なPython環境をゼロから構築します。すべて無料・商用利用可能なツールのみです。

## 前提

- macOS / Linux / Windows (WSL2)
- ターミナル（コマンドライン）が使えること

## micromambaのインストール

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 確認
micromamba --version
```

**ポイント:** micromambaはcondaの軽量版です。condaと同じパッケージが使えますが、インストール・解決が格段に速いです。

## Python環境の作成

```bash
# environment.yml を作成
cat << 'EOF' > environment.yml
name: crc-proteomics
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - numpy>=1.24
  - pandas>=2.0
  - scipy>=1.11
  - scikit-learn>=1.3
  - matplotlib>=3.7
  - seaborn>=0.13
  - statsmodels>=0.14
  - requests>=2.31
  - openpyxl>=3.1
  - pip
  - pip:
    - gseapy>=1.0
EOF

# 環境作成（5-10分）
micromamba create -f environment.yml -y

# 有効化
micromamba activate crc-proteomics
```

## 動作確認

```bash
python -c "
import numpy, pandas, scipy, sklearn, matplotlib, seaborn, statsmodels
print('All packages OK!')
print(f'NumPy: {numpy.__version__}')
print(f'Pandas: {pandas.__version__}')
print(f'SciPy: {scipy.__version__}')
"
```

## 注意点

- `micromamba activate` は毎回ターミナル起動時に実行が必要
- Windowsの場合はWSL2上で実行すること
- `CondaHTTPError` が出たらネットワーク接続を確認

## まとめ

micromamba + Python 3.11 で、DIA-MSプロテオミクス解析に必要な環境が構築できました。Perseus相当の統計解析がすべてPythonで実行できます。

> 前回: [#0 論文紹介](article-00-introduction.md)
> 次回: [#2 データ取得](article-02-data.md) — ProteomeXchangeからデータをダウンロード

#バイオインフォマティクス #プロテオミクス #Python #labcode
