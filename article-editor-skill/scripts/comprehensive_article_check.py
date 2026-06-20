#!/usr/bin/env python3
"""
記事の包括的品質チェック・book-design連携スクリプト
article-editorスキルのメイン実行スクリプト
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Tuple

# 同じディレクトリの他のスクリプトをインポート
sys.path.append(os.path.dirname(__file__))

try:
    from check_article_structure import generate_structure_report
    from check_figure_consistency import generate_figure_report
except ImportError as e:
    print(f"スクリプトインポートエラー: {e}")
    sys.exit(1)

class ArticleEditor:
    """記事エディタのメインクラス"""

    def __init__(self, blog_dir: str, paper_figures_dir: str = None):
        self.blog_dir = blog_dir
        self.paper_figures_dir = paper_figures_dir
        self.report = {}

    def run_comprehensive_check(self) -> Dict[str, Any]:
        """包括的な品質チェックを実行"""
        print("=== Article Editor: 包括的品質チェック開始 ===")

        # 1. 記事構成チェック
        print("\n[1/3] 記事構成の分析...")
        structure_report = generate_structure_report(self.blog_dir)

        # 2. 図表チェック
        print("\n[2/3] 図表整合性の分析...")
        figure_report = generate_figure_report(self.blog_dir, self.paper_figures_dir)

        # 3. 統合評価
        print("\n[3/3] 総合評価の実施...")
        comprehensive_report = self._generate_comprehensive_evaluation(
            structure_report, figure_report
        )

        self.report = comprehensive_report
        return comprehensive_report

    def _generate_comprehensive_evaluation(
        self, structure_report: Dict, figure_report: Dict
    ) -> Dict[str, Any]:
        """構成チェックと図チェックの結果を統合評価"""

        # エラーを統合
        all_errors = []
        all_warnings = []
        all_info = []

        # 構成レポートから問題を抽出
        if structure_report.get("status") == "success":
            issues = structure_report.get("issues", {})
            all_errors.extend(issues.get("errors", []))
            all_warnings.extend(issues.get("warnings", []))
            all_info.extend(issues.get("info", []))

        # 図レポートから問題を抽出
        if figure_report.get("status") == "success":
            issues = figure_report.get("issues", {})
            all_errors.extend(issues.get("errors", []))
            all_warnings.extend(issues.get("warnings", []))
            all_info.extend(issues.get("info", []))

        # 品質判定
        quality_assessment = self._assess_overall_quality(all_errors, all_warnings, all_info)

        # book-design連携判定
        book_design_action = self._determine_book_design_action(
            all_errors, all_warnings, quality_assessment
        )

        return {
            "timestamp": datetime.now().isoformat(),
            "status": "success",
            "structure_report": structure_report,
            "figure_report": figure_report,
            "consolidated_issues": {
                "errors": all_errors,
                "warnings": all_warnings,
                "info": all_info,
                "total_count": len(all_errors) + len(all_warnings) + len(all_info)
            },
            "quality_assessment": quality_assessment,
            "book_design_action": book_design_action,
            "recommendations": self._generate_action_recommendations(
                all_errors, all_warnings, quality_assessment
            )
        }

    def _assess_overall_quality(
        self, errors: List, warnings: List, info: List
    ) -> Dict[str, Any]:
        """全体的な品質を評価"""

        error_count = len(errors)
        warning_count = len(warnings)
        info_count = len(info)

        # 品質スコア算出（100点満点）
        score = 100
        score -= error_count * 20  # エラー1個で20点減点
        score -= warning_count * 5   # 警告1個で5点減点
        score -= info_count * 1      # 情報1個で1点減点
        score = max(0, score)

        # 品質レベル判定
        if score >= 90:
            quality_level = "excellent"
            quality_description = "出版品質: 優秀"
        elif score >= 75:
            quality_level = "good"
            quality_description = "出版品質: 良好"
        elif score >= 60:
            quality_level = "acceptable"
            quality_description = "出版品質: 許容範囲（要改善）"
        elif score >= 40:
            quality_level = "poor"
            quality_description = "出版品質: 不良（修正必要）"
        else:
            quality_level = "critical"
            quality_description = "出版品質: 重大問題（大幅修正必要）"

        return {
            "score": score,
            "level": quality_level,
            "description": quality_description,
            "error_count": error_count,
            "warning_count": warning_count,
            "info_count": info_count,
            "publishable": score >= 75 and error_count == 0
        }

    def _determine_book_design_action(
        self, errors: List, warnings: List, quality_assessment: Dict
    ) -> Dict[str, Any]:
        """book-designへのアクション（差し戻し判定）を決定"""

        # 差し戻し基準
        send_back_to_book_design = (
            quality_assessment["error_count"] > 0 or  # エラーが1個でもある
            quality_assessment["warning_count"] > 5 or  # 警告が5個超
            quality_assessment["score"] < 60  # スコアが60未満
        )

        if send_back_to_book_design:
            # 差し戻し時の詳細フィードバック生成
            feedback = self._generate_book_design_feedback(errors, warnings)
            action = "send_back"
            message = "重要な問題が見つかりました。book-designでの修正が必要です。"
        else:
            # 軽微な修正のみ
            feedback = self._generate_minor_fix_suggestions(warnings)
            action = "minor_fixes"
            message = "軽微な改善提案があります。article-editorで対応可能です。"

        return {
            "action": action,
            "message": message,
            "feedback": feedback,
            "send_back": send_back_to_book_design
        }

    def _generate_book_design_feedback(
        self, errors: List, warnings: List
    ) -> Dict[str, Any]:
        """book-designへの具体的フィードバックを生成"""

        feedback = {
            "critical_issues": [],
            "structural_improvements": [],
            "content_suggestions": [],
            "technical_fixes": []
        }

        # エラーを重要な問題として分類
        for error in errors:
            issue_type = error.get("type", "unknown")
            description = error.get("description", "")
            files = error.get("files", [])

            critical_issue = {
                "type": issue_type,
                "description": description,
                "affected_files": files,
                "priority": "high"
            }

            if issue_type in ["duplicate_title", "broken_internal_link"]:
                critical_issue["suggested_action"] = "記事構成の見直し・統合を検討"
                feedback["structural_improvements"].append(critical_issue)
            elif issue_type in ["missing_image", "figure_reference_mismatch"]:
                critical_issue["suggested_action"] = "図表の修正・追加"
                feedback["technical_fixes"].append(critical_issue)
            else:
                feedback["critical_issues"].append(critical_issue)

        # 警告を改善提案として分類
        for warning in warnings:
            issue_type = warning.get("type", "unknown")
            description = warning.get("description", "")

            if issue_type in ["numbering_gap", "similar_structure"]:
                feedback["structural_improvements"].append({
                    "type": issue_type,
                    "description": description,
                    "priority": "medium"
                })
            elif issue_type in ["terminology_inconsistency", "poor_alt_text"]:
                feedback["content_suggestions"].append({
                    "type": issue_type,
                    "description": description,
                    "priority": "medium"
                })

        return feedback

    def _generate_minor_fix_suggestions(self, warnings: List) -> Dict[str, Any]:
        """軽微な修正提案を生成"""

        suggestions = {
            "terminology_fixes": [],
            "alt_text_improvements": [],
            "formatting_adjustments": []
        }

        for warning in warnings:
            issue_type = warning.get("type", "")
            if issue_type == "terminology_inconsistency":
                suggestions["terminology_fixes"].append(warning)
            elif issue_type == "poor_alt_text":
                suggestions["alt_text_improvements"].append(warning)
            else:
                suggestions["formatting_adjustments"].append(warning)

        return suggestions

    def _generate_action_recommendations(
        self, errors: List, warnings: List, quality_assessment: Dict
    ) -> List[str]:
        """実行すべきアクションの推奨事項を生成"""

        recommendations = []

        if quality_assessment["publishable"]:
            recommendations.append("✅ 出版品質基準を満たしています")
        else:
            recommendations.append("❌ 出版前に修正が必要です")

        if errors:
            recommendations.append(f"🔴 重要: {len(errors)}件のエラーを修正してください")

        if warnings:
            recommendations.append(f"🟡 警告: {len(warnings)}件の改善提案があります")

        if quality_assessment["score"] < 75:
            recommendations.append("📚 book-designでの構成見直しを推奨")

        if quality_assessment["score"] >= 90:
            recommendations.append("🌟 優秀な品質です！publication-editorへ進めます")

        return recommendations

    def save_report(self, output_file: str = None):
        """レポートを保存"""
        if not output_file:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"article_editor_report_{timestamp}.json"

        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(self.report, f, ensure_ascii=False, indent=2)

        print(f"\n📄 レポートを保存しました: {output_file}")

    def print_summary(self):
        """結果のサマリーを表示"""
        if not self.report:
            print("レポートが生成されていません")
            return

        quality = self.report.get("quality_assessment", {})
        action = self.report.get("book_design_action", {})
        recommendations = self.report.get("recommendations", [])

        print("\n" + "="*60)
        print("📊 ARTICLE EDITOR - 品質チェック結果")
        print("="*60)

        print(f"\n🎯 品質スコア: {quality.get('score', 0)}/100")
        print(f"📈 品質レベル: {quality.get('description', 'N/A')}")
        print(f"📚 出版可能: {'はい' if quality.get('publishable') else 'いいえ'}")

        print(f"\n📋 問題数:")
        issues = self.report.get("consolidated_issues", {})
        print(f"  🔴 エラー: {len(issues.get('errors', []))}件")
        print(f"  🟡 警告: {len(issues.get('warnings', []))}件")
        print(f"  🔵 情報: {len(issues.get('info', []))}件")

        print(f"\n🎬 推奨アクション:")
        print(f"  📝 {action.get('message', 'N/A')}")
        if action.get('send_back'):
            print("  🔄 book-designに差し戻しが必要です")
        else:
            print("  ✨ article-editorで調整可能です")

        print(f"\n💡 推奨事項:")
        for rec in recommendations:
            print(f"  • {rec}")

def main():
    """メイン関数"""
    if len(sys.argv) < 2:
        print("使用法: python comprehensive_article_check.py <blog_directory> [paper_figures_directory] [output_file]")
        print("\n例:")
        print("  python comprehensive_article_check.py Article/blog/")
        print("  python comprehensive_article_check.py Article/blog/ data/paper_figures/ report.json")
        sys.exit(1)

    blog_directory = sys.argv[1]
    paper_figures_directory = sys.argv[2] if len(sys.argv) > 2 else None
    output_file = sys.argv[3] if len(sys.argv) > 3 else None

    # Article Editorを実行
    editor = ArticleEditor(blog_directory, paper_figures_directory)

    try:
        report = editor.run_comprehensive_check()
        editor.print_summary()
        editor.save_report(output_file)

        # 差し戻し判定に基づく終了ステータス
        action = report.get("book_design_action", {})
        if action.get("send_back"):
            print("\n🔄 book-designでの修正が必要です。")
            sys.exit(2)  # 差し戻し必要
        else:
            print("\n✅ 記事品質チェック完了。")
            sys.exit(0)  # 正常終了

    except Exception as e:
        print(f"\n❌ エラーが発生しました: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()