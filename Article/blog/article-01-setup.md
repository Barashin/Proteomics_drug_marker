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
    - adjustText>=0.8
    - pymzml>=2.5
EOF

# 環境作成（5-10分）
micromamba create -f environment.yml -y

# 有効化
micromamba activate crc-proteomics
```

## 動作確認

```bash
python -c "
import numpy as np, pandas as pd, matplotlib as plt, seaborn as sns
import scipy, sklearn, statsmodels, adjustText, pymzml, gseapy, openpyxl
print('✅ 基本パッケージ動作確認完了')
print(f'NumPy: {np.__version__}')
print(f'Pandas: {pd.__version__}')
print(f'matplotlib: {plt.__version__}')
print(f'seaborn: {sns.__version__}')
print(f'SciPy: {scipy.__version__}')
print(f'scikit-learn: {sklearn.__version__}')
print(f'adjustText: 利用可能')
print(f'pymzml: {pymzml.__version__}')
print(f'gseapy: {gseapy.__version__}')
"
```

**期待される出力:**
```
✅ 基本パッケージ動作確認完了
NumPy: 1.26.4
Pandas: 2.2.3
matplotlib: 3.10.9
seaborn: 0.13.2
SciPy: 1.12.0
scikit-learn: 1.8.0
adjustText: 利用可能
pymzml: 2.6
gseapy: 1.2.1
```

## 注意点

- `micromamba activate` は毎回ターミナル起動時に実行が必要
- Windowsの場合はWSL2上で実行すること
- `CondaHTTPError` が出たらネットワーク接続を確認

## GPU加速環境の構築【推奨】

より高速な解析のため、GPU加速環境を構築します。RTX 4070等のNVIDIA GPUが必要です。

### GPU環境設定ファイル作成

```bash
cat << 'EOF' > environment-gpu.yml
name: crc-proteomics-gpu
channels:
  - conda-forge
  - bioconda
  - pytorch
  - nvidia
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
  # GPU support - WSL2 compatible
  - pytorch=2.0.1
  - torchvision=0.15.2
  - pytorch-cuda=11.8
  # Proteomics tools
  - openms=3.1.0
  - pip
  - pip:
    - peptdeep>=1.4.0
    - pyprophet>=2.2.0
    - gseapy>=1.0
EOF
```

### GPU環境作成

```bash
# GPU環境作成（10-15分）
micromamba env create -f environment-gpu.yml -y

# 有効化
micromamba activate crc-proteomics-gpu

# PyTorchをpip経由で再インストール（WSL2互換性向上）
pip uninstall torch torchvision -y
pip install torch==2.0.1 torchvision==0.15.2 --index-url https://download.pytorch.org/whl/cu118

# プロテオミクスツール追加
pip install peptdeep>=1.4.0 pyprophet>=2.2.0
```

### GPU環境動作確認

```bash
python -c "
# OpenMS関連パッケージ確認
import pyopenms
print(f'✅ pyopenms: {pyopenms.__version__}')

import peptdeep
print(f'✅ peptdeep (AlphaPeptDeep): {peptdeep.__version__}')

import pyprophet
print('✅ pyprophet: 利用可能')

# PyTorch確認
import torch
print(f'✅ PyTorch: {torch.__version__}')
print(f'CUDA利用可能: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'GPU: {torch.cuda.get_device_name(0)}')
    memory_gb = torch.cuda.get_device_properties(0).total_memory / 1024**3
    print(f'GPU メモリ: {memory_gb:.1f} GB')
    print('🎉 GPU ready!')
else:
    print('⚠️ CPU版で動作（深層学習は可能だが低速）')
"
```

**期待される出力（GPU利用可能な場合）:**
```
✅ pyopenms: 3.5.0
✅ peptdeep (AlphaPeptDeep): 1.4.2
✅ pyprophet: 利用可能
✅ PyTorch: 2.11.0+cu118
CUDA利用可能: True
GPU: NVIDIA GeForce RTX 4070
GPU メモリ: 12.0 GB
🎉 GPU ready!
```

**CPU版の場合（深層学習は可能）:**
```
✅ pyopenms: 3.5.0
✅ peptdeep (AlphaPeptDeep): 1.4.2
✅ pyprophet: 利用可能
✅ PyTorch: 2.11.0+cpu
CUDA利用可能: False
⚠️ CPU版で動作（深層学習は可能だが低速）
```

### トラブルシューティング

**1. CUDA初期化エラー**
```bash
# WSL2でよくある問題。PyTorchをpip経由で再インストール
pip uninstall torch torchvision -y
pip install torch==2.0.1 torchvision==0.15.2 --index-url https://download.pytorch.org/whl/cu118
```

**2. 依存関係競合**
```bash
# 環境を削除して再作成
micromamba env remove -n crc-proteomics-gpu -y
micromamba env create -f environment-gpu.yml -y
```

## パフォーマンス比較（実測結果）

| 処理段階 | CPU環境 | GPU環境 | 加速倍率 | GPU使用率 |
|---------|---------|---------|----------|-----------|
| 電荷予測 | ~2分 | 21秒 | **6倍** | 95% |
| RT予測 | ~5分 | 21秒 | **15倍** | 90% |
| MS2予測 | ~10分 | 20秒 | **30倍** | 85% |
| **総時間** | **~20分** | **~2分** | **10倍高速** | - |

**RTX 4070環境での実測値。Toyota et al. 2025パラメータ（771,738プリカーサー）で検証済み。**

## 🎯 2つの環境の使い分け

| 環境 | 用途 | 検出性能 | 推奨場面 |
|------|------|----------|----------|
| **crc-proteomics** | **Sage基礎解析** | 2,110タンパク質 | 初学者、軽量解析 |
| **crc-proteomics-gpu** | **OpenMS深層学習** | 19,981タンパク質（9.5倍） | 高性能解析、論文級結果 |

## 📓 Jupyter Notebook の実行

### 基本解析（推奨・初心者向け）
```bash
# 基本環境を有効化
micromamba activate crc-proteomics

# Jupyter起動
jupyter notebook notebooks/

# ブラウザで http://localhost:8888 が開きます
```

### 深層学習解析（高性能・上級者向け）
```bash
# GPU環境を有効化
micromamba activate crc-proteomics-gpu

# Jupyter起動
jupyter notebook notebooks/

# notebook_10_openms.ipynb 以降で使用
```

### 💡 実行のコツ
- **notebook_00〜09**: `crc-proteomics` で実行
- **notebook_10〜16**: `crc-proteomics-gpu` で実行
- 環境切り替えは Jupyter を再起動してから

## 🚨 追加トラブルシューティング

### 3. PyTorchの __version__ エラー
```bash
# 破損したPyTorchを修正
micromamba activate crc-proteomics-gpu
pip uninstall torch torchvision -y
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
```

### 4. Jupyter Notebookが見つからない
```bash
# 基本環境にJupyterを追加インストール
micromamba activate crc-proteomics
pip install jupyter nbconvert
```

## まとめ

- **基本環境**: 2,110タンパク質解析（Sage）
- **GPU環境**: 19,981タンパク質解析（OpenMS深層学習）
- **使い分け**: Notebook番号で環境選択

micromamba + 2つの環境で、初心者から上級者まで対応できるDIA-MSプロテオミクス解析環境が構築できました。

> 前回: [#0 はじめに](article-00-introduction.md) — プロジェクト概要、全体像
> 次回: [#2a データ取得](article-02a-data-acquisition.md) — ProteomeXchangeからデータをダウンロード

#バイオインフォマティクス #プロテオミクス #Python #GPU #CUDA #labcode
