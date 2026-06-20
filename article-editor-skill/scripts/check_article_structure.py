#!/usr/bin/env python3
"""
記事全体の構成・整合性チェックスクリプト
article-editorスキルで使用される
"""

import os
import re
import json
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Any

def analyze_article_structure(blog_dir: str) -> Dict[str, Any]:
    """記事全体の構造を解析"""
    blog_path = Path(blog_dir)
    if not blog_path.exists():
        return {"error": f"ディレクトリが見つかりません: {blog_dir}"}

    articles = []
    issues = []

    # 記事ファイルを収集
    article_files = sorted([f for f in blog_path.glob("article-*.md")])

    for article_file in article_files:
        try:
            with open(article_file, 'r', encoding='utf-8') as f:
                content = f.read()

            # 基本情報を抽出
            title_match = re.search(r'^#\s+(.+)$', content, re.MULTILINE)
            title = title_match.group(1).strip() if title_match else "タイトルなし"

            # セクション情報を抽出
            sections = re.findall(r'^##\s+(.+)$', content, re.MULTILINE)

            # 図表情報を抽出
            figures = re.findall(r'!\[([^\]]*)\]\(([^)]+)\)', content)
            tables = len(re.findall(r'^\|', content, re.MULTILINE))

            # コードブロック情報
            code_blocks = len(re.findall(r'```', content)) // 2

            # 参照情報
            internal_links = re.findall(r'\[([^\]]+)\]\(article-([^)]+)\)', content)

            articles.append({
                "filename": article_file.name,
                "title": title,
                "sections": sections,
                "figures": figures,
                "table_count": tables,
                "code_block_count": code_blocks,
                "internal_links": internal_links,
                "word_count": len(content.split()),
                "content": content
            })

        except Exception as e:
            issues.append(f"ファイル読み込みエラー: {article_file.name} - {str(e)}")

    return {
        "articles": articles,
        "issues": issues,
        "total_count": len(articles)
    }

def check_consistency_issues(articles: List[Dict]) -> List[Dict[str, Any]]:
    """構成の一貫性をチェック"""
    issues = []

    # 1. ファイル命名の一貫性チェック
    article_numbers = []
    for article in articles:
        filename = article["filename"]
        # article-XX または article-XXa のようなパターンを抽出
        match = re.match(r'article-(\d+)([a-z]?)-', filename)
        if match:
            number = int(match.group(1))
            sub = match.group(2) or ""
            article_numbers.append((number, sub, filename))
        else:
            issues.append({
                "type": "naming_inconsistency",
                "severity": "warning",
                "description": f"命名規則に合わないファイル: {filename}",
                "files": [filename]
            })

    # 連続性チェック
    article_numbers.sort()
    for i, (num, sub, filename) in enumerate(article_numbers[1:], 1):
        prev_num, prev_sub, prev_filename = article_numbers[i-1]
        if num - prev_num > 1 and not sub and not prev_sub:
            issues.append({
                "type": "numbering_gap",
                "severity": "warning",
                "description": f"記事番号に飛びがあります: {prev_filename} → {filename}",
                "files": [prev_filename, filename]
            })

    # 2. 内容の重複チェック
    title_counts = Counter([article["title"] for article in articles])
    for title, count in title_counts.items():
        if count > 1:
            duplicate_files = [art["filename"] for art in articles if art["title"] == title]
            issues.append({
                "type": "duplicate_title",
                "severity": "error",
                "description": f"重複するタイトル: '{title}'",
                "files": duplicate_files
            })

    # 3. セクション構造の一貫性チェック
    section_patterns = defaultdict(list)
    for article in articles:
        pattern = " → ".join(article["sections"])
        section_patterns[pattern].append(article["filename"])

    # 似たようなセクション構成の記事を特定
    similar_patterns = [(pattern, files) for pattern, files in section_patterns.items() if len(files) > 1]
    if similar_patterns:
        issues.append({
            "type": "similar_structure",
            "severity": "info",
            "description": "類似の構成を持つ記事グループが見つかりました",
            "details": similar_patterns
        })

    # 4. 内部リンクの整合性チェック
    all_filenames = {article["filename"].replace(".md", "") for article in articles}
    for article in articles:
        for link_text, link_target in article["internal_links"]:
            target_file = f"article-{link_target}.md"
            if target_file.replace(".md", "") not in all_filenames:
                issues.append({
                    "type": "broken_internal_link",
                    "severity": "error",
                    "description": f"存在しない記事への内部リンク: '{link_text}' → {target_file}",
                    "files": [article["filename"]]
                })

    # 5. 記事長のバランスチェック
    word_counts = [article["word_count"] for article in articles]
    avg_words = sum(word_counts) / len(word_counts)
    for article in articles:
        if article["word_count"] < avg_words * 0.3:
            issues.append({
                "type": "short_article",
                "severity": "warning",
                "description": f"極端に短い記事: {article['filename']} ({article['word_count']}語)",
                "files": [article["filename"]]
            })
        elif article["word_count"] > avg_words * 3:
            issues.append({
                "type": "long_article",
                "severity": "warning",
                "description": f"極端に長い記事: {article['filename']} ({article['word_count']}語)",
                "files": [article["filename"]]
            })

    return issues

def check_technical_consistency(articles: List[Dict]) -> List[Dict[str, Any]]:
    """技術的内容の整合性をチェック"""
    issues = []

    # 1. 用語の統一チェック
    terminology_variants = [
        (["proteomic", "proteomics"], "プロテオミクス関連用語"),
        (["DIA", "dia", "DIA-MS"], "DIA用語"),
        (["OpenMS", "openms", "OpenSwath"], "OpenMS関連用語"),
        (["Sage", "sage", "SAGE"], "Sage用語"),
    ]

    for article in articles:
        content = article["content"].lower()
        for variants, term_type in terminology_variants:
            found_variants = [var for var in variants if var.lower() in content]
            if len(found_variants) > 1:
                issues.append({
                    "type": "terminology_inconsistency",
                    "severity": "warning",
                    "description": f"用語の表記ゆれ ({term_type}): {found_variants}",
                    "files": [article["filename"]]
                })

    # 2. 図表参照の整合性
    for article in articles:
        content = article["content"]
        # 図の参照を検索
        fig_refs = re.findall(r'図\s*(\d+)', content)
        fig_numbers = [int(ref) for ref in fig_refs]

        # 連続性チェック
        if fig_numbers:
            fig_numbers.sort()
            for i, num in enumerate(fig_numbers[1:], 1):
                if num - fig_numbers[i-1] > 1:
                    issues.append({
                        "type": "figure_numbering_gap",
                        "severity": "warning",
                        "description": f"図番号に飛びがあります: {article['filename']}",
                        "files": [article["filename"]]
                    })

    return issues

def generate_structure_report(blog_dir: str) -> Dict[str, Any]:
    """記事構造の完全な分析レポートを生成"""
    print("記事構造を解析中...")
    structure_data = analyze_article_structure(blog_dir)

    if "error" in structure_data:
        return structure_data

    articles = structure_data["articles"]

    print("整合性をチェック中...")
    consistency_issues = check_consistency_issues(articles)

    print("技術的整合性をチェック中...")
    technical_issues = check_technical_consistency(articles)

    all_issues = consistency_issues + technical_issues

    # 重要度別に問題を分類
    errors = [issue for issue in all_issues if issue["severity"] == "error"]
    warnings = [issue for issue in all_issues if issue["severity"] == "warning"]
    info = [issue for issue in all_issues if issue["severity"] == "info"]

    # 統計情報
    stats = {
        "total_articles": len(articles),
        "total_words": sum(article["word_count"] for article in articles),
        "total_figures": sum(len(article["figures"]) for article in articles),
        "total_tables": sum(article["table_count"] for article in articles),
        "total_code_blocks": sum(article["code_block_count"] for article in articles),
    }

    return {
        "status": "success",
        "statistics": stats,
        "articles": articles,
        "issues": {
            "errors": errors,
            "warnings": warnings,
            "info": info,
            "total_count": len(all_issues)
        },
        "recommendation": generate_recommendations(errors, warnings, stats)
    }

def generate_recommendations(errors: List, warnings: List, stats: Dict) -> Dict[str, Any]:
    """問題解決の推奨事項を生成"""
    recommendations = {
        "critical_actions": [],
        "improvements": [],
        "book_design_feedback": []
    }

    if errors:
        recommendations["critical_actions"].append(
            "エラーレベルの問題があります。book-designスキルでの修正が必要です。"
        )
        for error in errors:
            recommendations["book_design_feedback"].append({
                "type": error["type"],
                "description": error["description"],
                "affected_files": error.get("files", [])
            })

    if len(warnings) > 5:
        recommendations["improvements"].append(
            "多数の警告があります。記事構成の見直しを推奨します。"
        )

    if stats["total_articles"] > 30:
        recommendations["improvements"].append(
            f"記事数が{stats['total_articles']}個と多いです。統合できる記事がないか検討してください。"
        )

    return recommendations

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("使用法: python check_article_structure.py <blog_directory>")
        sys.exit(1)

    blog_directory = sys.argv[1]
    report = generate_structure_report(blog_directory)

    print(json.dumps(report, ensure_ascii=False, indent=2))