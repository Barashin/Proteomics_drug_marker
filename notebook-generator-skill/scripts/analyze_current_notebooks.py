#!/usr/bin/env python3
"""
現在のnotebook状況分析スクリプト
notebook-generatorスキルで使用
"""

import os
import json
import re
from pathlib import Path
from typing import Dict, List, Any

def analyze_articles_structure(blog_dir: str) -> Dict[str, Any]:
    """記事構造を分析"""
    blog_path = Path(blog_dir)
    articles = []

    article_files = sorted(blog_path.glob("article-*.md"))

    for article_file in article_files:
        with open(article_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # コードブロックを抽出
        code_blocks = re.findall(r'```(\w+)?\n(.*?)```', content, re.DOTALL)
        python_blocks = [block[1] for block in code_blocks if block[0] in ['python', 'py', '']]
        bash_blocks = [block[1] for block in code_blocks if block[0] in ['bash', 'shell']]

        # 記事番号を抽出
        filename = article_file.name
        match = re.match(r'article-(\d+)([a-z]?)-', filename)
        if match:
            number = int(match.group(1))
            sub = match.group(2)
        else:
            number = 999
            sub = ""

        articles.append({
            "filename": filename,
            "number": number,
            "sub": sub,
            "python_code_blocks": len(python_blocks),
            "bash_code_blocks": len(bash_blocks),
            "total_code_lines": sum(len(block.split('\n')) for block in python_blocks + bash_blocks),
            "has_code": len(python_blocks) + len(bash_blocks) > 0
        })

    return {
        "articles": articles,
        "total_articles": len(articles),
        "articles_with_code": sum(1 for a in articles if a["has_code"])
    }

def analyze_notebooks_structure(notebooks_dir: str) -> Dict[str, Any]:
    """既存notebook構造を分析"""
    notebooks_path = Path(notebooks_dir)
    notebooks = []

    if not notebooks_path.exists():
        return {"notebooks": [], "total_notebooks": 0}

    notebook_files = sorted(notebooks_path.glob("*.ipynb"))

    for notebook_file in notebook_files:
        try:
            with open(notebook_file, 'r', encoding='utf-8') as f:
                nb_data = json.load(f)

            # セル数を分析
            cells = nb_data.get("cells", [])
            markdown_cells = sum(1 for cell in cells if cell.get("cell_type") == "markdown")
            code_cells = sum(1 for cell in cells if cell.get("cell_type") == "code")

            # ステップ番号を抽出
            filename = notebook_file.name
            match = re.match(r'step_(\d+)', filename)
            if match:
                step_number = int(match.group(1))
            else:
                step_number = 999

            notebooks.append({
                "filename": filename,
                "step_number": step_number,
                "markdown_cells": markdown_cells,
                "code_cells": code_cells,
                "total_cells": len(cells)
            })

        except Exception as e:
            notebooks.append({
                "filename": filename,
                "error": str(e)
            })

    return {
        "notebooks": notebooks,
        "total_notebooks": len(notebooks)
    }

def suggest_article_notebook_mapping(articles: List[Dict], notebooks: List[Dict]) -> Dict[str, Any]:
    """記事とnotebookの対応マッピングを提案"""

    # 記事をグループ化（同じ番号の記事をまとめる）
    article_groups = {}
    for article in articles:
        number = article["number"]
        if number not in article_groups:
            article_groups[number] = []
        article_groups[number].append(article)

    # 提案マッピング
    suggested_mapping = {}
    missing_notebooks = []

    for number in sorted(article_groups.keys()):
        group = article_groups[number]
        notebook_name = f"notebook-{number:02d}.ipynb"

        # コードを含む記事があるかチェック
        has_code = any(article["has_code"] for article in group)

        if has_code:
            suggested_mapping[notebook_name] = [article["filename"] for article in group]

            # 既存のnotebookがあるかチェック
            existing = any(nb.get("step_number") == number for nb in notebooks)
            if not existing:
                missing_notebooks.append(notebook_name)

    return {
        "suggested_mapping": suggested_mapping,
        "missing_notebooks": missing_notebooks,
        "total_notebooks_needed": len(suggested_mapping)
    }

def generate_analysis_report(blog_dir: str, notebooks_dir: str) -> Dict[str, Any]:
    """包括的な分析レポートを生成"""

    print("記事構造を分析中...")
    articles_analysis = analyze_articles_structure(blog_dir)

    print("既存notebook構造を分析中...")
    notebooks_analysis = analyze_notebooks_structure(notebooks_dir)

    print("記事とnotebookの対応関係を分析中...")
    mapping_suggestion = suggest_article_notebook_mapping(
        articles_analysis["articles"],
        notebooks_analysis["notebooks"]
    )

    # 統計情報
    stats = {
        "articles_total": articles_analysis["total_articles"],
        "articles_with_code": articles_analysis["articles_with_code"],
        "notebooks_current": notebooks_analysis["total_notebooks"],
        "notebooks_needed": mapping_suggestion["total_notebooks_needed"],
        "notebooks_missing": len(mapping_suggestion["missing_notebooks"])
    }

    return {
        "timestamp": "2026-05-01",
        "statistics": stats,
        "articles_analysis": articles_analysis,
        "notebooks_analysis": notebooks_analysis,
        "mapping_suggestion": mapping_suggestion,
        "recommendations": generate_recommendations(stats, mapping_suggestion)
    }

def generate_recommendations(stats: Dict, mapping: Dict) -> List[str]:
    """推奨事項を生成"""
    recommendations = []

    if stats["notebooks_missing"] > 0:
        recommendations.append(f"🔴 {stats['notebooks_missing']}個のnotebookが不足しています")

    if stats["notebooks_current"] < stats["notebooks_needed"]:
        recommendations.append("📝 記事とnotebookの対応関係を整理する必要があります")

    if stats["articles_with_code"] < stats["articles_total"] * 0.5:
        recommendations.append("⚠️ コードを含まない記事が多いです（notebookは不要）")

    recommendations.append("✅ チェーン実行のための依存関係設定を推奨")
    recommendations.append("📋 全notebookでの共通ライブラリ統一を推奨")

    return recommendations

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("使用法: python analyze_current_notebooks.py <blog_dir> <notebooks_dir>")
        sys.exit(1)

    blog_directory = sys.argv[1]
    notebooks_directory = sys.argv[2]

    report = generate_analysis_report(blog_directory, notebooks_directory)

    print("\n" + "="*60)
    print("📊 NOTEBOOK 現状分析レポート")
    print("="*60)

    stats = report["statistics"]
    print(f"\n📈 統計:")
    print(f"  📄 記事数: {stats['articles_total']}個")
    print(f"  🔢 コード含む記事: {stats['articles_with_code']}個")
    print(f"  📓 現在のnotebook: {stats['notebooks_current']}個")
    print(f"  🎯 必要なnotebook: {stats['notebooks_needed']}個")
    print(f"  ❌ 不足notebook: {stats['notebooks_missing']}個")

    print(f"\n💡 推奨事項:")
    for rec in report["recommendations"]:
        print(f"  • {rec}")

    print(f"\n📋 提案マッピング:")
    for notebook, articles in report["mapping_suggestion"]["suggested_mapping"].items():
        print(f"  📓 {notebook}")
        for article in articles:
            print(f"    ← {article}")

    # 詳細レポートをJSONで保存
    output_file = "notebook_analysis_report.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    print(f"\n📄 詳細レポート保存: {output_file}")