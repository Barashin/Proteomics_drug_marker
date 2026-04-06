#!/bin/bash
# run_all.sh - Toyota et al. 2025 再現解析パイプライン
# 全ステップを順番に実行する
# 番号はブログ記事 (article-XX) と対応
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=========================================="
echo "Toyota et al. 2025 再現解析パイプライン"
echo "=========================================="
echo ""

echo "=== Step 02: データ取得（article-02） ==="
python "$SCRIPT_DIR/step_02_download_data.py"
echo ""

echo "=== Step 03: mzMLファイル検査（article-03 前半） ==="
python "$SCRIPT_DIR/step_03_inspect_mzml.py"
echo ""

echo "=== Step 04a: sage によるDIA解析（article-03） ==="
python "$SCRIPT_DIR/step_04_dia_analysis_sage.py"
echo ""

echo "=== Step 04b: タンパク質マトリクス構築（article-03） ==="
python "$SCRIPT_DIR/step_04_build_protein_matrix.py"
echo ""

echo "=== Step 04c: データ前処理（article-04） ==="
python "$SCRIPT_DIR/step_04_preprocess.py"
echo ""

echo "=== Step 05: 全体像の可視化（article-05） ==="
python "$SCRIPT_DIR/step_05_overview_visualization.py"
echo ""

echo "=== Step 06: 差分発現解析（article-06） ==="
python "$SCRIPT_DIR/step_06_differential_expression.py"
echo ""

echo "=== Step 07: COSMIC解析（article-07） ==="
python "$SCRIPT_DIR/step_07_cosmic_analysis.py"
echo ""

echo "=== Step 08: ステージ別解析（article-08） ==="
python "$SCRIPT_DIR/step_08_stage_analysis.py"
echo ""

echo "=========================================="
echo "全ステップ完了！"
echo "結果は results/ ディレクトリに保存されています。"
echo "=========================================="
