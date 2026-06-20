#!/usr/bin/env python3
"""
不足notebookの自動生成スクリプト
チェーン実行可能な形式で生成
"""

import os
import re
import json
import nbformat as nbf
from pathlib import Path
from typing import List, Dict, Any

# 不足notebook一覧（現状分析結果から）
MISSING_NOTEBOOKS = [
    {
        "notebook": "notebook-00.ipynb",
        "articles": ["article-00-introduction.md"],
        "description": "DIA-MSプロテオミクス プロジェクト概要"
    },
    {
        "notebook": "notebook-01.ipynb",
        "articles": ["article-01-setup.md"],
        "description": "環境構築とセットアップ"
    },
    {
        "notebook": "notebook-02.ipynb",
        "articles": ["article-02a-data-acquisition.md", "article-02b-data-formats.md"],
        "description": "データ取得とファイル形式理解"
    },
    {
        "notebook": "notebook-03.ipynb",
        "articles": ["article-03a-convert-basics.md", "article-03b-convert-execution.md"],
        "description": "RAW→mzML変換"
    },
    {
        "notebook": "notebook-06.ipynb",
        "articles": ["article-06a-preprocess-basics.md", "article-06b-preprocess-imputation.md"],
        "description": "データ前処理と欠損値補完"
    },
    {
        "notebook": "notebook-09.ipynb",
        "articles": ["article-09-cosmic.md"],
        "description": "COSMICデータベースとの照合"
    }
]

def extract_code_from_article(article_path: str) -> List[Dict[str, str]]:
    """記事からコードブロックとMarkdownを抽出"""

    with open(article_path, 'r', encoding='utf-8') as f:
        content = f.read()

    sections = []
    current_section = ""

    # YAMLフロントマターをスキップ
    if content.startswith('---'):
        parts = content.split('---', 2)
        if len(parts) >= 3:
            content = parts[2]

    lines = content.split('\n')
    code_block = ""
    in_code_block = False
    code_type = ""

    for line in lines:
        if line.startswith('```'):
            if in_code_block:
                # コードブロック終了
                if code_block.strip():
                    sections.append({
                        "type": "code",
                        "language": code_type,
                        "content": code_block.strip()
                    })
                code_block = ""
                in_code_block = False
                code_type = ""
            else:
                # コードブロック開始
                if current_section.strip():
                    sections.append({
                        "type": "markdown",
                        "content": current_section.strip()
                    })
                    current_section = ""
                in_code_block = True
                code_type = line[3:].strip() if len(line) > 3 else "text"
        elif in_code_block:
            code_block += line + "\n"
        else:
            current_section += line + "\n"

    # 残りのセクションを追加
    if current_section.strip():
        sections.append({
            "type": "markdown",
            "content": current_section.strip()
        })

    return sections

def create_notebook_header(notebook_info: Dict, articles_dir: str) -> List[Dict]:
    """Notebookヘッダーセルを作成"""

    article_links = []
    for article in notebook_info["articles"]:
        article_path = f"../blog/{article}"
        title = article.replace(".md", "").replace("article-", "#").replace("-", " ")
        article_links.append(f"- [{title}]({article_path})")

    header_content = f"""# {notebook_info['description']}

> **対応するブログ記事:**
{chr(10).join(article_links)}
>
> このNotebookは上記ブログ記事のコードをセルごとに実行できるインタラクティブ版です。

---
"""

    return [{
        "type": "markdown",
        "content": header_content
    }]

def create_common_imports() -> Dict:
    """共通ライブラリのimportセル"""

    import_code = """# === 共通ライブラリのインポート ===
import os, glob, re, subprocess, sys, time, warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# 警告を抑制、matplotlib設定
warnings.filterwarnings('ignore')
plt.style.use('default')
%matplotlib inline

# プロジェクト共通パス設定
PROJECT_ROOT = Path("../").resolve()
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

# ディレクトリの確認・作成
for dir_path in [DATA_DIR, RESULTS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

print("✅ 環境設定完了")
print(f"📁 プロジェクトルート: {PROJECT_ROOT}")
print(f"📁 データディレクトリ: {DATA_DIR}")
print(f"📁 結果ディレクトリ: {RESULTS_DIR}")"""

    return {
        "type": "code",
        "content": import_code
    }

def create_dependency_check(notebook_info: Dict) -> Dict:
    """前ステップ依存関係チェックセル"""

    notebook_num = int(re.search(r'notebook-(\d+)', notebook_info["notebook"]).group(1))

    if notebook_num == 0:
        # 最初のnotebook
        dependency_code = '''print("📚 プロジェクト開始 - 依存関係チェックなし")'''
    else:
        prev_notebook = f"notebook-{notebook_num-1:02d}.ipynb"
        dependency_code = f'''# === 前ステップの依存関係チェック ===
prev_results = [
    RESULTS_DIR / "step_{notebook_num-1:02d}_completed.flag",
    # 必要に応じて他のファイルも追加
]

missing_files = [f for f in prev_results if not f.exists()]

if missing_files:
    print("⚠️ 前のステップの実行が必要です:")
    for f in missing_files:
        print(f"  ❌ 不足: {{f}}")
    print(f"\\n👉 先に {prev_notebook} を実行してください")
    raise FileNotFoundError("前ステップの出力ファイルが見つかりません")
else:
    print("✅ 依存関係チェック完了 - 前ステップの出力確認済み")'''

    return {
        "type": "code",
        "content": dependency_code
    }

def create_completion_flag(notebook_info: Dict) -> Dict:
    """完了フラグ作成セル"""

    notebook_num = int(re.search(r'notebook-(\d+)', notebook_info["notebook"]).group(1))

    completion_code = f'''# === ステップ完了フラグの作成 ===
completion_flag = RESULTS_DIR / "step_{notebook_num:02d}_completed.flag"
completion_flag.touch()

print("✅ {notebook_info['description']} 完了")
print(f"🚩 完了フラグ作成: {{completion_flag}}")
print(f"👉 次のステップ: notebook-{notebook_num+1:02d}.ipynb")'''

    return {
        "type": "code",
        "content": completion_code
    }

def generate_notebook(notebook_info: Dict, articles_dir: str, output_dir: str):
    """単一notebookの生成"""

    nb = nbf.v4.new_notebook()
    cells = []

    # 1. ヘッダーセル
    header_cells = create_notebook_header(notebook_info, articles_dir)
    for cell_info in header_cells:
        if cell_info["type"] == "markdown":
            cells.append(nbf.v4.new_markdown_cell(cell_info["content"]))

    # 2. 共通ライブラリimport
    import_cell = create_common_imports()
    cells.append(nbf.v4.new_code_cell(import_cell["content"]))

    # 3. 依存関係チェック
    dependency_cell = create_dependency_check(notebook_info)
    cells.append(nbf.v4.new_code_cell(dependency_cell["content"]))

    # 4. 記事内容を結合
    all_sections = []
    for article in notebook_info["articles"]:
        article_path = Path(articles_dir) / article
        if article_path.exists():
            print(f"  処理中: {article}")
            sections = extract_code_from_article(str(article_path))

            # 記事区切りセクション
            if len(notebook_info["articles"]) > 1:
                separator = f"## {article.replace('.md', '').replace('article-', '#').replace('-', ' ')}"
                all_sections.append({
                    "type": "markdown",
                    "content": f"{separator}\\n\\n---\\n*このセクションは [{article}](../blog/{article}) の内容です*\\n"
                })

            all_sections.extend(sections)
        else:
            print(f"  ⚠️ ファイルが見つかりません: {article}")

    # セクションをセルに変換
    for section in all_sections:
        if section["type"] == "markdown":
            cells.append(nbf.v4.new_markdown_cell(section["content"]))
        elif section["type"] == "code":
            # Python以外のコードは適切に変換
            code_content = section["content"]
            if section.get("language") == "bash":
                code_content = f"# Bash コマンド\\n!{code_content}"
            elif section.get("language") == "yaml":
                # YAMLは文字列として出力するコードに変換
                code_content = f'# YAML設定\\nyaml_content = """\\n{code_content}\\n"""\\nprint("YAML設定:")\\nprint(yaml_content)'

            cells.append(nbf.v4.new_code_cell(code_content))

    # 5. 完了フラグセル
    completion_cell = create_completion_flag(notebook_info)
    cells.append(nbf.v4.new_code_cell(completion_cell["content"]))

    nb.cells = cells

    # notebook保存
    output_path = Path(output_dir) / notebook_info["notebook"]
    with open(output_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)

    print(f"✅ 生成完了: {output_path}")

def main():
    import sys

    if len(sys.argv) != 3:
        print("使用法: python generate_missing_notebooks.py <articles_dir> <output_dir>")
        sys.exit(1)

    articles_dir = sys.argv[1]
    output_dir = sys.argv[2]

    # 出力ディレクトリの作成
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print("🚀 不足notebookの生成開始")
    print(f"📁 記事ディレクトリ: {articles_dir}")
    print(f"📁 出力ディレクトリ: {output_dir}")
    print(f"📊 生成対象: {len(MISSING_NOTEBOOKS)}個")

    # 各notebookを生成
    for i, notebook_info in enumerate(MISSING_NOTEBOOKS, 1):
        print(f"\\n[{i}/{len(MISSING_NOTEBOOKS)}] {notebook_info['notebook']}")
        try:
            generate_notebook(notebook_info, articles_dir, output_dir)
        except Exception as e:
            print(f"❌ エラー: {e}")

    print(f"\\n🎉 生成完了: {len(MISSING_NOTEBOOKS)}個のnotebook")
    print("\\n📋 チェーン実行方法:")
    print("  1. notebook-00.ipynb から順番に実行")
    print("  2. 各notebookは前ステップの完了を確認してから実行")
    print("  3. エラーが出た場合は前のnotebookを確認")

if __name__ == "__main__":
    main()