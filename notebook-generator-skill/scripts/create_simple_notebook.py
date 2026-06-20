#!/usr/bin/env python3
"""
シンプルなnotebook生成（nbformat不要版）
チェーン実行対応
"""

import os
import re
import json
from pathlib import Path
from typing import List, Dict

def create_notebook_json(cells: List[Dict]) -> Dict:
    """Jupyter notebook JSONを生成"""
    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "name": "python",
                "version": "3.11.0"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 4
    }

def create_markdown_cell(content: str) -> Dict:
    """Markdownセルを作成"""
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": content.split('\n')
    }

def create_code_cell(content: str) -> Dict:
    """Codeセルを作成"""
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": content.split('\n')
    }

def extract_code_blocks(article_path: str) -> List[str]:
    """記事からPythonコードブロックを抽出"""
    with open(article_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # YAMLフロントマターをスキップ
    if content.startswith('---'):
        parts = content.split('---', 2)
        if len(parts) >= 3:
            content = parts[2]

    # ```python または ```py または ``` のコードブロックを抽出
    code_blocks = []
    pattern = r'```(?:python|py)?\n(.*?)```'
    matches = re.findall(pattern, content, re.DOTALL)

    for match in matches:
        cleaned = match.strip()
        if cleaned and not cleaned.startswith('bash') and not cleaned.startswith('yaml'):
            code_blocks.append(cleaned)

    return code_blocks

def create_chain_notebook(notebook_num: int, description: str, articles: List[str], articles_dir: str) -> Dict:
    """チェーン実行対応notebookを作成"""

    cells = []

    # 1. ヘッダーセル
    header = f"""# Notebook {notebook_num:02d}: {description}

> **対応記事:** {', '.join(articles)}
>
> このNotebookは記事のコードをインタラクティブに実行できます。

---
"""
    cells.append(create_markdown_cell(header))

    # 2. 共通ライブラリ
    common_imports = '''# 共通ライブラリのインポート
import os, glob, re, subprocess, sys, time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
plt.style.use('default')

# プロジェクト設定
PROJECT_ROOT = Path("..").resolve()
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
for d in [DATA_DIR, RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

print("✅ 環境設定完了")'''
    cells.append(create_code_cell(common_imports))

    # 3. 依存関係チェック
    if notebook_num > 0:
        dependency_check = f'''# 前ステップの依存関係チェック
prev_flag = RESULTS_DIR / "step_{notebook_num-1:02d}_completed.flag"
if not prev_flag.exists():
    print(f"⚠️ 前のステップを先に実行してください: notebook-{notebook_num-1:02d}.ipynb")
    raise FileNotFoundError("前ステップの完了フラグが見つかりません")
else:
    print("✅ 依存関係チェック完了")'''
    else:
        dependency_check = '''print("📚 プロジェクト開始 - 依存関係なし")'''

    cells.append(create_code_cell(dependency_check))

    # 4. 記事のコード
    for article in articles:
        article_path = Path(articles_dir) / article
        if article_path.exists():
            # 記事セクションヘッダー
            section_header = f"## {article}"
            cells.append(create_markdown_cell(section_header))

            # コードブロック抽出
            code_blocks = extract_code_blocks(str(article_path))
            for i, code in enumerate(code_blocks):
                cells.append(create_code_cell(code))

    # 5. 完了フラグ
    completion_flag = f'''# ステップ完了
completion_flag = RESULTS_DIR / "step_{notebook_num:02d}_completed.flag"
completion_flag.touch()
print(f"✅ Step {notebook_num:02d} 完了: {description}")
print(f"👉 次のステップ: notebook-{notebook_num+1:02d}.ipynb")'''
    cells.append(create_code_cell(completion_flag))

    return create_notebook_json(cells)

# 不足notebook情報
MISSING_NOTEBOOKS = [
    {"num": 0, "desc": "プロジェクト概要", "articles": ["article-00-introduction.md"]},
    {"num": 1, "desc": "環境構築", "articles": ["article-01-setup.md"]},
    {"num": 2, "desc": "データ取得", "articles": ["article-02a-data-acquisition.md", "article-02b-data-formats.md"]},
    {"num": 3, "desc": "データ変換", "articles": ["article-03a-convert-basics.md", "article-03b-convert-execution.md"]},
    {"num": 6, "desc": "前処理", "articles": ["article-06a-preprocess-basics.md", "article-06b-preprocess-imputation.md"]},
    {"num": 9, "desc": "COSMIC照合", "articles": ["article-09-cosmic.md"]}
]

def main():
    import sys

    if len(sys.argv) != 3:
        print("使用法: python create_simple_notebook.py <articles_dir> <output_dir>")
        sys.exit(1)

    articles_dir = sys.argv[1]
    output_dir = sys.argv[2]

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print("🚀 チェーン実行対応notebook生成開始")

    for info in MISSING_NOTEBOOKS:
        notebook_name = f"notebook-{info['num']:02d}.ipynb"
        print(f"  📓 {notebook_name}: {info['desc']}")

        # notebook JSON作成
        nb_json = create_chain_notebook(
            info['num'],
            info['desc'],
            info['articles'],
            articles_dir
        )

        # ファイル保存
        output_path = Path(output_dir) / notebook_name
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(nb_json, f, indent=2, ensure_ascii=False)

        print(f"    ✅ 生成完了: {output_path}")

    print(f"\n🎉 {len(MISSING_NOTEBOOKS)}個のnotebook生成完了")
    print("\n📋 チェーン実行方法:")
    print("  1. notebook-00.ipynb から順番に実行")
    print("  2. 各セルを上から順に実行")
    print("  3. 依存関係エラーが出たら前のnotebookを確認")

if __name__ == "__main__":
    main()