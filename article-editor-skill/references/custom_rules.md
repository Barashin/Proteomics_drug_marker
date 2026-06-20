# カスタム検証ルール

プロテオミクス・創薬マーカー記事プロジェクト専用の品質チェックルール

## 用語統一ルール

### 必須統一用語

| 正式表記 | 避けるべき表記 | 備考 |
|---------|---------------|------|
| DIA-MS | dia-ms, DIA MS, dia | Data Independent Acquisition Mass Spectrometry |
| OpenMS | openms, OpenSwath, openswath | オープンソースMS解析ツール |
| Sage | sage, SAGE | 高速ペプチド検索エンジン |
| プロテオミクス | proteomics, プロテオーム解析 | 日本語記事では統一 |
| 質量分析 | MS, Mass Spec | 日本語記事では正式名称 |
| ペプチド | peptide, ペプタイド | 日本語表記統一 |
| タンパク質 | protein, プロテイン | 日本語記事では統一 |

### ファイル名・パス表記ルール

| 項目 | ルール | 例 |
|------|-------|-----|
| データファイル | スネークケース | `protein_matrix.csv` |
| 図ファイル | 記事番号付き | `fig_05_pca_plot.png` |
| スクリプト | step番号付き | `step_05_openms.py` |
| 結果ディレクトリ | 小文字スラッシュ区切り | `results/figures/` |

## 記事構成ルール

### 記事番号規則

- **基本**: `article-XX-title.md` （XXは2桁番号）
- **サブ記事**: `article-XXa-title.md`, `article-XXb-title.md`
- **連続性**: 番号は連続的、大きな飛びは避ける

### セクション構造規則

各記事は以下の構造を基本とする：

```markdown
# タイトル

## はじめに
（背景・目的の簡潔な説明）

## 方法/実装
（具体的な手順・コード）

## 結果
（出力・図表）

## まとめ
（要点・次ステップへの橋渡し）
```

### 必須要素

- **前記事への言及**: 流れの明確化
- **次記事への予告**: シリーズの一貫性
- **図表番号**: 連続的な番号付け
- **コード例**: 実行可能な形で提供

## 図表品質ルール

### 図のファイル命名

```
fig_[記事番号]_[内容]_[種類].png
例: fig_05_pca_plot.png, fig_08_volcano_plot.png
```

### 必須Alt Text要素

- 図の種類（PCA plot, Volcano plot など）
- データの概要（サンプル数、条件など）
- 主要な知見

**良い例**: `PCA plot showing separation of 3 cancer stages (n=120 samples)`
**悪い例**: `plot`, `figure1`, `graph`

### キャプション要素

1. **図の説明**（何を示すか）
2. **データ詳細**（サンプル数、条件）
3. **解釈**（主要な発見）
4. **技術詳細**（必要に応じて）

## コード品質ルール

### インポート統一

```python
# 必須インポート順序
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# プロジェクト固有
import openms
from sage import search
```

### パス参照統一

```python
# データパス
DATA_DIR = "data/"
RESULTS_DIR = "results/"
FIGURES_DIR = "results/figures/"

# 相対パス禁止項目
× "../data/file.csv"
○ "data/file.csv"
```

### 関数命名規則

```python
# 解析関数
def analyze_protein_expression(df):
def generate_pca_plot(data, output_path):
def calculate_differential_proteins(df1, df2):

# ユーティリティ関数
def load_protein_matrix(file_path):
def save_figure(fig, filename):
```

## 品質基準

### エラーレベル（即修正必要）

- 存在しないファイル参照
- 重複タイトル
- 用語統一違反（必須統一用語）
- 図番号の重複・欠番
- 実行不可能なコード例

### 警告レベル（改善推奨）

- 記事番号の飛び（1つまでOK）
- Alt textが15文字未満
- セクション構造の不統一
- 図キャプションの不足
- 前後記事への言及不足

### 情報レベル（最適化提案）

- 類似内容記事の統合可能性
- 図の解像度・見やすさ改善
- コードコメントの充実
- 参考文献の追加提案

## チェック自動化

### 必須チェック項目

```bash
# 用語統一チェック
grep -r "dia-ms\|DIA MS" Article/blog/

# ファイル存在チェック
find Article/blog/ -name "*.md" -exec grep -l "!\[.*\](" {} \;

# 図番号チェック
grep -r "図[0-9]" Article/blog/
```

### カスタムバリデーション

プロジェクト固有の検証スクリプト `scripts/validate_proteomics_rules.py` でルール違反を検出。

## 除外設定

### チェック除外項目

- **引用文**: 論文からの直接引用は原文保持
- **固有名詞**: ツール名・データベース名は原則維持
- **コードコメント**: 英語コメントは統一対象外
- **URL・パス**: 外部リンクは変更禁止