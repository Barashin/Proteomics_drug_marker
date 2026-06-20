#!/usr/bin/env python3
"""
図表の整合性・視覚的比較スクリプト
article-editorスキルで使用される
"""

import os
import re
import json
from pathlib import Path
from typing import Dict, List, Tuple, Any
import hashlib

def extract_figure_info(blog_dir: str) -> Dict[str, Any]:
    """記事から図表情報を抽出"""
    blog_path = Path(blog_dir)
    figure_data = {}
    issues = []

    article_files = sorted(blog_path.glob("article-*.md"))

    for article_file in article_files:
        try:
            with open(article_file, 'r', encoding='utf-8') as f:
                content = f.read()

            # 図の情報を抽出
            figures = re.findall(r'!\[([^\]]*)\]\(([^)]+)\)', content)

            # 図の参照を抽出
            fig_references = re.findall(r'図\s*(\d+)', content)

            # 図のキャプション情報
            fig_captions = []
            for alt_text, image_path in figures:
                caption_context = extract_figure_caption(content, image_path)
                fig_captions.append({
                    "alt_text": alt_text,
                    "image_path": image_path,
                    "caption_context": caption_context
                })

            figure_data[article_file.name] = {
                "figures": figures,
                "fig_references": fig_references,
                "fig_captions": fig_captions,
                "total_figures": len(figures)
            }

        except Exception as e:
            issues.append(f"図情報抽出エラー: {article_file.name} - {str(e)}")

    return {
        "figure_data": figure_data,
        "issues": issues
    }

def extract_figure_caption(content: str, image_path: str) -> str:
    """画像周辺のキャプション情報を抽出"""
    # 画像の前後のテキストを取得
    image_pattern = rf'!\[[^\]]*\]\({re.escape(image_path)}\)'

    match = re.search(image_pattern, content)
    if not match:
        return ""

    start_pos = match.start()
    end_pos = match.end()

    # 前後100文字を取得
    before_text = content[max(0, start_pos-100):start_pos].strip()
    after_text = content[end_pos:end_pos+100].strip()

    return f"前: {before_text} | 後: {after_text}"

def check_figure_consistency(figure_data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """図の整合性をチェック"""
    issues = []

    # 1. 図ファイルの存在確認
    for article_name, data in figure_data.items():
        for alt_text, image_path in data["figures"]:
            # 相対パス解決
            if not image_path.startswith('http'):
                full_path = Path(image_path)
                if not full_path.exists():
                    issues.append({
                        "type": "missing_image",
                        "severity": "error",
                        "description": f"画像ファイルが見つかりません: {image_path}",
                        "article": article_name,
                        "image_path": image_path
                    })

    # 2. 図番号の整合性チェック
    for article_name, data in figure_data.items():
        fig_refs = [int(ref) for ref in data["fig_references"] if ref.isdigit()]
        actual_figures = data["total_figures"]

        if fig_refs:
            max_ref = max(fig_refs)
            if max_ref > actual_figures:
                issues.append({
                    "type": "figure_reference_mismatch",
                    "severity": "warning",
                    "description": f"図の参照番号({max_ref})が実際の図数({actual_figures})を超えています",
                    "article": article_name
                })

    # 3. 重複する画像ファイルのチェック
    image_paths = []
    for article_name, data in figure_data.items():
        for alt_text, image_path in data["figures"]:
            image_paths.append((image_path, article_name))

    path_counts = {}
    for path, article in image_paths:
        if path not in path_counts:
            path_counts[path] = []
        path_counts[path].append(article)

    for path, articles in path_counts.items():
        if len(articles) > 1:
            issues.append({
                "type": "duplicate_image",
                "severity": "info",
                "description": f"同じ画像が複数の記事で使用されています: {path}",
                "articles": articles
            })

    # 4. Alt textの品質チェック
    for article_name, data in figure_data.items():
        for alt_text, image_path in data["figures"]:
            if not alt_text or len(alt_text.strip()) < 3:
                issues.append({
                    "type": "poor_alt_text",
                    "severity": "warning",
                    "description": f"図のalt textが不十分です: {image_path}",
                    "article": article_name,
                    "image_path": image_path
                })

    return issues

def compare_with_paper_figures(figure_data: Dict[str, Any], paper_figures_dir: str = None) -> List[Dict[str, Any]]:
    """論文の元図と比較（論文図が利用可能な場合）"""
    issues = []

    if not paper_figures_dir or not Path(paper_figures_dir).exists():
        issues.append({
            "type": "no_paper_reference",
            "severity": "info",
            "description": "論文の元図との比較ができません（参考図ディレクトリが見つかりません）",
        })
        return issues

    # 論文図ファイルを取得
    paper_figures = list(Path(paper_figures_dir).glob("*.png")) + \
                   list(Path(paper_figures_dir).glob("*.jpg")) + \
                   list(Path(paper_figures_dir).glob("*.jpeg"))

    if not paper_figures:
        issues.append({
            "type": "no_paper_figures",
            "severity": "info",
            "description": f"論文図ディレクトリに図が見つかりません: {paper_figures_dir}",
        })
        return issues

    # 記事の図と論文図の類似性をチェック（ファイル名ベース）
    article_images = []
    for article_name, data in figure_data.items():
        for alt_text, image_path in data["figures"]:
            if not image_path.startswith('http'):
                article_images.append((Path(image_path).name, article_name, image_path))

    paper_image_names = [fig.name for fig in paper_figures]

    for img_name, article_name, img_path in article_images:
        # 部分的な名前の一致をチェック
        img_stem = Path(img_name).stem.lower()
        similar_papers = [p for p in paper_image_names if img_stem in p.lower() or p.lower() in img_stem]

        if similar_papers:
            issues.append({
                "type": "potential_paper_match",
                "severity": "info",
                "description": f"論文図との類似性を確認してください: {img_name}",
                "article": article_name,
                "image_path": img_path,
                "similar_paper_figures": similar_papers
            })

    return issues

def generate_figure_report(blog_dir: str, paper_figures_dir: str = None) -> Dict[str, Any]:
    """図の包括的な分析レポートを生成"""
    print("図情報を抽出中...")
    figure_info = extract_figure_info(blog_dir)

    if figure_info["issues"]:
        return {
            "status": "error",
            "errors": figure_info["issues"]
        }

    figure_data = figure_info["figure_data"]

    print("図の整合性をチェック中...")
    consistency_issues = check_figure_consistency(figure_data)

    print("論文図との比較を実行中...")
    paper_comparison_issues = compare_with_paper_figures(figure_data, paper_figures_dir)

    all_issues = consistency_issues + paper_comparison_issues

    # 重要度別に分類
    errors = [issue for issue in all_issues if issue["severity"] == "error"]
    warnings = [issue for issue in all_issues if issue["severity"] == "warning"]
    info = [issue for issue in all_issues if issue["severity"] == "info"]

    # 統計情報
    total_figures = sum(data["total_figures"] for data in figure_data.values())
    articles_with_figures = len([data for data in figure_data.values() if data["total_figures"] > 0])

    stats = {
        "total_figures": total_figures,
        "articles_with_figures": articles_with_figures,
        "articles_without_figures": len(figure_data) - articles_with_figures,
        "avg_figures_per_article": total_figures / len(figure_data) if figure_data else 0
    }

    return {
        "status": "success",
        "statistics": stats,
        "figure_data": figure_data,
        "issues": {
            "errors": errors,
            "warnings": warnings,
            "info": info,
            "total_count": len(all_issues)
        },
        "recommendations": generate_figure_recommendations(errors, warnings, stats)
    }

def generate_figure_recommendations(errors: List, warnings: List, stats: Dict) -> Dict[str, Any]:
    """図に関する推奨事項を生成"""
    recommendations = {
        "critical_actions": [],
        "improvements": [],
        "visual_checks": []
    }

    if errors:
        recommendations["critical_actions"].append(
            "図に関する重要な問題があります。修正が必要です。"
        )

    if stats["articles_without_figures"] > stats["articles_with_figures"] * 0.5:
        recommendations["improvements"].append(
            "図のない記事が多いです。視覚的な説明の追加を検討してください。"
        )

    if stats["avg_figures_per_article"] > 10:
        recommendations["improvements"].append(
            "記事あたりの図が多いです。必要な図に絞ることを検討してください。"
        )

    # 視覚的チェックの推奨
    potential_matches = [issue for issue in warnings if issue.get("type") == "potential_paper_match"]
    if potential_matches:
        recommendations["visual_checks"].append(
            "論文図との類似性がある図が見つかりました。視覚的に比較してください。"
        )

    return recommendations

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("使用法: python check_figure_consistency.py <blog_directory> [paper_figures_directory]")
        sys.exit(1)

    blog_directory = sys.argv[1]
    paper_figures_directory = sys.argv[2] if len(sys.argv) > 2 else None

    report = generate_figure_report(blog_directory, paper_figures_directory)
    print(json.dumps(report, ensure_ascii=False, indent=2))