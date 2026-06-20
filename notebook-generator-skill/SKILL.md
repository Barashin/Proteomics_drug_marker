---
name: notebook-generator
description: 記事→Notebook自動生成スキル - blog記事からJupyter Notebookを生成し、記事とnotebookの対応関係を整理する。記事のMarkdown部分は説明セルに、コードブロックは実行可能セルに変換。記事番号とnotebook番号を統一し、チェーン実行可能な形で依存関係を管理。図はinline表示、ライブラリ統一、前ステップ出力の自動読み込みを実現。labcode-publisherの記事→notebook化ステップで使用。「notebook作成」「記事をnotebookに」「notebook-generator」「記事対応」等でトリガー。
---

# Notebook Generator - 記事→Notebook自動生成

blog記事からJupyter Notebookを自動生成し、実行可能な学習環境を構築する。

## 機能概要

### 1. 記事→Notebook変換
- **Markdownセル**: 記事の説明部分をそのまま変換
- **Codeセル**: コードブロックを実行可能形式で抽出
- **図の処理**: inline表示でnotebook内に埋め込み

### 2. 統合記事対応（新方針）
- **統合記事 → Notebook**: article-02-data.md → notebook-02.ipynb
- **セクション統合**: 複数記事の内容を自然なセクション構成に
- **Zenn出版対応**: チャプター構成を考慮した見出し調整
- **依存関係の自動処理**: 前ステップ出力の読み込みコード追加

### 3. 実行環境の統一
- **共通ライブラリ**: 全notebookで同一のimport文
- **パス設定**: 相対パスで統一、プロジェクト構造に対応
- **エラーハンドリング**: ファイル存在確認、gracefulな処理

## 実行手順

### Step 1: プロジェクト構造の分析

```bash
python scripts/analyze_project_structure.py Article/blog/ Article/notebooks/
```

**分析項目**:
- 既存記事の番号体系とコード内容
- 現在のnotebook状況
- 記事間の依存関係（内部リンク、ファイル参照）
- 必要なライブラリとインポート

### Step 2: 記事→Notebook変換マッピング

記事の細分化状況を分析してnotebook統合方針を決定：

```python
# 記事とnotebook名完全対応マッピング（新方針）
conversion_map = {
    "notebook-00-introduction.ipynb": ["article-00-introduction.md"],
    "notebook-01-setup.ipynb": ["article-01-setup.md"],
    "notebook-02-data.ipynb": ["article-02-data.md"],  # 統合記事対応
    "notebook-03-convert.ipynb": ["article-03-convert.md"],  # 統合記事対応
    "notebook-04-sage.ipynb": ["article-04-sage.md"],  # 統合記事対応
    "notebook-05-data-exploration.ipynb": ["article-05-data-exploration.md"],  # 統合記事対応
    "notebook-06-preprocess.ipynb": ["article-06-preprocess.md"],  # 統合記事対応
    "notebook-07-visualization.ipynb": ["article-07-visualization.md"],  # 統合記事対応
    "notebook-08-differential.ipynb": ["article-08-differential.md"],  # 統合記事対応
    "notebook-09-cosmic.ipynb": ["article-09-cosmic.md"],
    "notebook-10-openms.ipynb": ["article-10-openms.md"],
    "notebook-11-openms-stage.ipynb": ["article-11-openms-stage.md"],  # 統合記事対応
    "notebook-12-openms-visualization.ipynb": ["article-12-openms-visualization.md"],  # 統合記事対応
    "notebook-13-openms-differential.ipynb": ["article-13-openms-differential.md"],  # 統合記事対応
    "notebook-14-openms-cosmic.ipynb": ["article-14-openms-cosmic.md"],
    "notebook-15-openms-stage.ipynb": ["article-15-openms-stage.md"],  # 統合記事対応
    "notebook-16-comparison.ipynb": ["article-16-comparison.md"],  # 統合記事対応
    "notebook-17-conclusion.ipynb": ["article-17-conclusion.md"]
}
```

### Step 3: Notebook生成

各notebookを順次生成：

```bash
python scripts/generate_notebooks.py --input Article/blog/ --output Article/notebooks/ --config conversion_config.json
```

**生成内容**:
- ヘッダーセル（対応記事へのリンク、概要）
- ライブラリimportセル（共通化）
- 記事内容のmarkdownセル
- コードブロックのcodeセル
- 結果確認セル（図表示、統計出力）

### Step 4: 依存関係の設定

前ステップ出力を読み込むコードを自動挿入：

```python
# 例：notebook-05で前ステップの結果を読み込む
# 自動生成されるコード
if os.path.exists("../results/protein_matrix.csv"):
    df = pd.read_csv("../results/protein_matrix.csv", index_col=0)
    print(f"前ステップの結果を読み込み: {df.shape}")
else:
    print("前ステップを先に実行してください: notebook-04.ipynb")
    raise FileNotFoundError("Required input file not found")
```

### Step 5: 実行テスト

生成されたnotebook群の実行テスト：

```bash
python scripts/test_notebook_execution.py Article/notebooks/
```

## 変換ルール

### Markdownセル変換

**記事のセクション → Markdownセル**:
```markdown
# 記事のタイトル → # Markdownセル
## セクション見出し → ## Markdownセル
説明文 → そのまま転記
> 注釈ブロック → > 注釈ブロック
```

**対応記事リンクの挿入**:
```markdown
> **対応するブログ記事**: [#4 Sage-proteomicsでDIA解析](../blog/article-04a-sage-fundamentals.md)
>
> このNotebookはブログ記事 #4 のコードをセルごとに実行できるインタラクティブ版です。
```

### Codeセル変換

**コードブロック抽出ルール**:
```python
# 記事の```python ブロック → Jupyter Codeセル
# 記事の```bash ブロック → !bash コマンドとしてCodeセル
# 記事の```yaml ブロック → ファイル出力用Codeセル
```

**共通ライブラリのimport**:
```python
# 全notebookの冒頭に挿入される標準import
import os, glob, re, subprocess, sys, time, warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
warnings.filterwarnings('ignore')
plt.style.use('default')
%matplotlib inline
```

### パス設定の統一

**相対パス基準の設定**:
```python
# プロジェクトルートからの相対パス設定
PROJECT_ROOT = Path("../").resolve()
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

# ディレクトリ存在確認
for dir_path in [DATA_DIR, RESULTS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)
```

## 高度な機能

### 図表の自動処理

**inline表示の実装**:
```python
# 図の生成・表示コード例
plt.figure(figsize=(10, 6))
# ... plotting code ...
plt.title("タンパク質発現量の分布")
plt.show()  # inline表示

# ファイル保存も併用
plt.savefig(RESULTS_DIR / "figures" / "protein_distribution.png", dpi=300, bbox_inches='tight')
```

### エラーハンドリング

**Gracefulな処理**:
```python
try:
    # 必要なファイルの確認
    required_files = [
        DATA_DIR / "raw" / "protein_matrix.csv",
        DATA_DIR / "reference" / "cosmic_genes.txt"
    ]

    for file_path in required_files:
        if not file_path.exists():
            print(f"⚠️  ファイルが見つかりません: {file_path}")
            print("前のnotebookを実行してください")
            raise FileNotFoundError(f"Required file missing: {file_path}")

    # メイン処理
    process_data()

except Exception as e:
    print(f"❌ エラーが発生しました: {e}")
    print("詳細なエラー情報:")
    import traceback
    traceback.print_exc()
```

### 統合記事からのNotebook生成（新方針）

**統合記事 → Notebook + Zenn対応**:
```python
# article-04-sage.md → notebook-04.ipynb + Zenn chapter構成
def create_unified_notebook(unified_article, zenn_structure):
    """統合記事から実行可能なNotebook生成"""

    # 統合記事の各セクションをNotebookセルに変換
    sections = parse_unified_article_sections(unified_article)

    notebook_cells = []
    for section in sections:
        # Zenn章構成を考慮した見出しレベル調整
        if section['type'] == 'heading':
            level = adjust_heading_for_zenn(section['level'])
            notebook_cells.append(create_markdown_cell(f"{'#' * level} {section['title']}"))
        elif section['type'] == 'code':
            notebook_cells.append(create_code_cell(section['content']))

    return create_notebook_json(notebook_cells)
```

## 品質管理

### 自動テスト

**Notebook実行テスト**:
```python
def test_notebook_execution(notebook_path):
    """Notebookが最後まで実行できるかテスト"""
    import nbformat
    from nbconvert.preprocessors import ExecutePreprocessor

    with open(notebook_path) as f:
        nb = nbformat.read(f, as_version=4)

    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
    try:
        ep.preprocess(nb, {'metadata': {'path': str(notebook_path.parent)}})
        return True, "OK"
    except Exception as e:
        return False, str(e)
```

**依存関係チェック**:
```python
def check_dependencies(notebooks_dir):
    """前ステップの出力ファイルが存在するかチェック"""
    for notebook in sorted(notebooks_dir.glob("notebook-*.ipynb")):
        # notebook内の依存ファイルをチェック
        pass
```

### 記事との同期

**差分チェック機能**:
```bash
python scripts/check_article_notebook_sync.py Article/blog/ Article/notebooks/
```

記事が更新された場合、対応するnotebookも自動更新。

## 使用例

### 基本的な変換

```bash
# 全記事をnotebookに変換
python scripts/generate_notebooks.py \
  --input Article/blog/ \
  --output Article/notebooks/ \
  --mode full

# 特定記事のみ変換
python scripts/generate_notebooks.py \
  --input Article/blog/article-04*.md \
  --output Article/notebooks/ \
  --mode single --target notebook-04
```

### 実行テスト

```bash
# 全notebook実行テスト
python scripts/test_all_notebooks.py Article/notebooks/

# 順次実行（依存関係考慮）
python scripts/run_notebooks_chain.py Article/notebooks/
```

### 記事更新時の同期

```bash
# 記事更新検出→notebook自動更新
python scripts/sync_articles_notebooks.py \
  --watch Article/blog/ \
  --notebooks Article/notebooks/
```

## 他スキルとの連携

### labcode-publisher ワークフロー統合

**Step A7.5: notebook-generator**（book-design後に挿入）
- book-designで記事が確定
- notebook-generatorで実行可能版作成
- article-editorで記事とnotebook両方の品質チェック

### code-debug 連携

生成されたnotebookの実行エラーを検出・修正：
- notebook実行時のエラーログ収集
- code-debugによる自動修正
- 修正版notebookの再生成

### code-explanator 連携

初学者向けnotebook版の生成：
- コード解説をmarkdownセルに追加
- 段階的実行のガイド挿入
- エラー対処法の説明

Notebook Generatorは技術記事を実践的な学習教材に変換し、読者が実際にコードを実行しながら学べる環境を提供する。