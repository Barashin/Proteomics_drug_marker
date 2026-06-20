#!/usr/bin/env python3
"""
記事文字数チェックスクリプト
article-editorスキルで使用される
"""

import os
import re
import json
from pathlib import Path
from typing import Dict, List, Any

# 記事タイプ別文字数基準
LENGTH_GUIDELINES = {
    "introduction": {"min": 1500, "max": 3000, "type": "導入記事"},
    "setup": {"min": 1500, "max": 3000, "type": "環境構築"},
    "implementation": {"min": 2000, "max": 5000, "type": "実装記事"},
    "analysis": {"min": 2500, "max": 4000, "type": "解析記事"},
    "comparison": {"min": 3000, "max": 6000, "type": "比較記事"},
    "conclusion": {"min": 2000, "max": 4000, "type": "まとめ記事"},
    "default": {"min": 2000, "max": 5000, "type": "標準記事"}
}

def classify_article_type(filename: str, title: str, content: str) -> str:
    """記事のタイプを分類"""

    filename_lower = filename.lower()
    title_lower = title.lower()
    content_lower = content.lower()

    # ファイル名による分類
    if "introduction" in filename_lower or "00" in filename_lower:
        return "introduction"
    elif "setup" in filename_lower or "01" in filename_lower:
        return "setup"
    elif "conclusion" in filename_lower or "17" in filename_lower:
        return "conclusion"
    elif "comparison" in filename_lower or "16" in filename_lower:
        return "comparison"

    # タイトル・内容による分類
    if any(word in title_lower for word in ["はじめに", "概要", "導入"]):
        return "introduction"
    elif any(word in title_lower for word in ["環境", "セットアップ", "構築"]):
        return "setup"
    elif any(word in title_lower for word in ["実装", "コード", "sage", "openms"]):
        return "implementation"
    elif any(word in title_lower for word in ["解析", "可視化", "統計", "差分"]):
        return "analysis"
    elif any(word in title_lower for word in ["比較", "評価", "性能"]):
        return "comparison"
    elif any(word in title_lower for word in ["まとめ", "結論", "総括"]):
        return "conclusion"

    return "default"

def analyze_article_length(article_path: str) -> Dict[str, Any]:
    """記事の詳細文字数分析"""

    with open(article_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # YAMLフロントマターを除去
    original_content = content
    if content.startswith('---'):
        parts = content.split('---', 2)
        if len(parts) >= 3:
            content = parts[2]

    # タイトル抽出
    title_match = re.search(r'^#\s+(.+)$', content, re.MULTILINE)
    title = title_match.group(1).strip() if title_match else "タイトルなし"

    # 総文字数
    total_chars = len(content)

    # コードブロックを抽出・除去
    code_blocks = re.findall(r'```.*?\n(.*?)```', content, re.DOTALL)
    code_chars = sum(len(block) for block in code_blocks)
    content_without_code = re.sub(r'```.*?```', '', content, flags=re.DOTALL)

    # 本文文字数（コード除く）
    text_chars = len(content_without_code)

    # 日本語文字数
    japanese_chars = len(re.findall(r'[ぁ-んァ-ヶー一-龯]', content_without_code))

    # 英数字文字数
    ascii_chars = len(re.findall(r'[a-zA-Z0-9]', content_without_code))

    # 図表数
    figures = len(re.findall(r'!\[([^\]]*)\]\(([^)]+)\)', content))
    tables = len(re.findall(r'^\|', content, re.MULTILINE))

    # 記事タイプを分類
    article_type = classify_article_type(Path(article_path).name, title, content)

    return {
        "filename": Path(article_path).name,
        "title": title,
        "article_type": article_type,
        "total_chars": total_chars,
        "text_chars": text_chars,  # コード除く本文
        "code_chars": code_chars,
        "japanese_chars": japanese_chars,
        "ascii_chars": ascii_chars,
        "figures": figures,
        "tables": tables,
        "code_ratio": code_chars / total_chars if total_chars > 0 else 0
    }

def check_length_compliance(article_info: Dict[str, Any]) -> List[Dict[str, Any]]:
    """文字数基準への準拠チェック"""

    issues = []
    article_type = article_info["article_type"]
    text_chars = article_info["text_chars"]

    # 基準を取得
    if article_type in LENGTH_GUIDELINES:
        guideline = LENGTH_GUIDELINES[article_type]
    else:
        guideline = LENGTH_GUIDELINES["default"]

    min_chars = guideline["min"]
    max_chars = guideline["max"]
    type_name = guideline["type"]

    # 文字数チェック
    if text_chars < 1000:
        issues.append({
            "type": "article_too_short",
            "severity": "error",
            "description": f"記事が短すぎます: {text_chars}文字（最低1,000文字必要）",
            "article": article_info["filename"],
            "expected": "1,000文字以上",
            "actual": f"{text_chars}文字"
        })
    elif text_chars < min_chars:
        issues.append({
            "type": "article_short",
            "severity": "warning",
            "description": f"{type_name}として短いです: {text_chars}文字（推奨{min_chars}-{max_chars}文字）",
            "article": article_info["filename"],
            "expected": f"{min_chars}-{max_chars}文字",
            "actual": f"{text_chars}文字"
        })
    elif text_chars > max_chars:
        if text_chars > 8000:
            issues.append({
                "type": "article_too_long",
                "severity": "warning",
                "description": f"記事が長すぎます: {text_chars}文字（分割検討推奨）",
                "article": article_info["filename"],
                "expected": "8,000文字以下",
                "actual": f"{text_chars}文字"
            })
        else:
            issues.append({
                "type": "article_long",
                "severity": "info",
                "description": f"{type_name}として長いです: {text_chars}文字（推奨{min_chars}-{max_chars}文字）",
                "article": article_info["filename"],
                "expected": f"{min_chars}-{max_chars}文字",
                "actual": f"{text_chars}文字"
            })

    # コード比率チェック
    code_ratio = article_info["code_ratio"]
    if code_ratio > 0.9:
        issues.append({
            "type": "too_much_code",
            "severity": "warning",
            "description": f"コード比率が高すぎます: {code_ratio:.1%}（説明不足の可能性）",
            "article": article_info["filename"],
            "expected": "90%以下",
            "actual": f"{code_ratio:.1%}"
        })
    elif code_ratio < 0.1 and article_info["article_type"] in ["implementation", "analysis"]:
        issues.append({
            "type": "too_little_code",
            "severity": "info",
            "description": f"コード比率が低いです: {code_ratio:.1%}（実装記事として）",
            "article": article_info["filename"],
            "expected": "10%以上",
            "actual": f"{code_ratio:.1%}"
        })

    return issues

def generate_length_report(blog_dir: str) -> Dict[str, Any]:
    """記事文字数の包括的レポート生成"""

    blog_path = Path(blog_dir)
    if not blog_path.exists():
        return {"error": f"ディレクトリが見つかりません: {blog_dir}"}

    print("記事文字数を分析中...")

    articles = []
    all_issues = []

    # 各記事を分析
    article_files = sorted(blog_path.glob("article-*.md"))
    for article_file in article_files:
        try:
            article_info = analyze_article_length(str(article_file))
            articles.append(article_info)

            # 文字数基準チェック
            issues = check_length_compliance(article_info)
            all_issues.extend(issues)

        except Exception as e:
            all_issues.append({
                "type": "analysis_error",
                "severity": "error",
                "description": f"記事分析エラー: {str(e)}",
                "article": article_file.name
            })

    # 統計情報
    if articles:
        total_text_chars = sum(a["text_chars"] for a in articles)
        avg_text_chars = total_text_chars // len(articles)

        text_lengths = [a["text_chars"] for a in articles]
        text_lengths.sort()
        median_text = text_lengths[len(text_lengths)//2]
        min_text = min(text_lengths)
        max_text = max(text_lengths)

        # タイプ別統計
        type_stats = {}
        for article in articles:
            art_type = article["article_type"]
            if art_type not in type_stats:
                type_stats[art_type] = []
            type_stats[art_type].append(article["text_chars"])

        for art_type, lengths in type_stats.items():
            type_stats[art_type] = {
                "count": len(lengths),
                "avg_length": sum(lengths) // len(lengths),
                "min_length": min(lengths),
                "max_length": max(lengths)
            }

        stats = {
            "total_articles": len(articles),
            "avg_text_chars": avg_text_chars,
            "median_text_chars": median_text,
            "min_text_chars": min_text,
            "max_text_chars": max_text,
            "type_distribution": type_stats
        }
    else:
        stats = {"total_articles": 0}

    # 重要度別問題分類
    errors = [issue for issue in all_issues if issue["severity"] == "error"]
    warnings = [issue for issue in all_issues if issue["severity"] == "warning"]
    info = [issue for issue in all_issues if issue["severity"] == "info"]

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
        "recommendations": generate_length_recommendations(stats, errors, warnings)
    }

def generate_length_recommendations(stats: Dict, errors: List, warnings: List) -> List[str]:
    """文字数に関する推奨事項生成"""

    recommendations = []

    if errors:
        recommendations.append(f"🔴 {len(errors)}件の重要な文字数問題があります")

    if warnings:
        recommendations.append(f"🟡 {len(warnings)}件の文字数調整推奨があります")

    if stats.get("total_articles", 0) > 0:
        avg_chars = stats.get("avg_text_chars", 0)
        if avg_chars < 2000:
            recommendations.append("📝 全体的に記事が短いです。内容の充実を検討してください")
        elif avg_chars > 6000:
            recommendations.append("📄 全体的に記事が長いです。分割を検討してください")

    recommendations.append("✅ 記事タイプに応じた文字数ガイドラインの遵守を推奨")

    return recommendations

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("使用法: python check_article_length.py <blog_directory>")
        sys.exit(1)

    blog_directory = sys.argv[1]
    report = generate_length_report(blog_directory)

    if "error" in report:
        print(f"エラー: {report['error']}")
        sys.exit(1)

    print("\n" + "="*60)
    print("📏 記事文字数チェック結果")
    print("="*60)

    stats = report["statistics"]
    if stats["total_articles"] > 0:
        print(f"\n📊 統計情報:")
        print(f"  記事数: {stats['total_articles']}個")
        print(f"  平均文字数: {stats['avg_text_chars']}文字")
        print(f"  中央値: {stats['median_text_chars']}文字")
        print(f"  範囲: {stats['min_text_chars']}-{stats['max_text_chars']}文字")

    issues = report["issues"]
    print(f"\n🔍 問題数:")
    print(f"  🔴 エラー: {len(issues['errors'])}件")
    print(f"  🟡 警告: {len(issues['warnings'])}件")
    print(f"  🔵 情報: {len(issues['info'])}件")

    if issues["errors"]:
        print(f"\n🔴 エラー詳細:")
        for error in issues["errors"][:5]:  # 最初の5件
            print(f"  • {error['description']}")

    if issues["warnings"]:
        print(f"\n🟡 警告詳細:")
        for warning in issues["warnings"][:5]:  # 最初の5件
            print(f"  • {warning['description']}")

    print(f"\n💡 推奨事項:")
    for rec in report["recommendations"]:
        print(f"  • {rec}")

    print(f"\n📄 詳細レポート保存: article_length_report.json")
    with open("article_length_report.json", 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)