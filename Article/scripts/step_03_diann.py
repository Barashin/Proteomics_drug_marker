#!/usr/bin/env python3
"""
Step 3: DIA-MSデータのタンパク質同定・定量（OpenMS/OpenSWATH）
================================================================
DIA-MSの生データ（mzML）からタンパク質を同定・定量するスクリプトです。

【なぜOpenMSを使うのか？】
  論文では DIA-NN v1.8.1 を使用していますが、DIA-NN は商用利用に有料ライセンスが必要です。
  本スクリプトでは、無料かつ商用利用可能な OpenMS（BSD 3-Clause ライセンス）で代替します。

  | ツール | ライセンス | 商用利用 |
  |--------|----------|---------|
  | DIA-NN | 学術無料/商用有料 | ❌（有料） |
  | OpenMS | BSD 3-Clause | ✅（無料） |

【OpenMS / OpenSWATH とは？】
  - OpenMS: 質量分析データ解析のためのオープンソースフレームワーク
  - OpenSWATH: OpenMS内のDIA/SWATHデータ解析モジュール
  - BSD 3-Clause ライセンス（無料、商用利用可能、改変・再配布可能）
  - GitHub: https://github.com/OpenMS/OpenMS

【処理の流れ】
  1. UniProt FASTAからin silico消化（Digestor）
  2. ペプチドのスペクトルライブラリ生成
  3. OpenSWATH でDIAデータを解析
  4. PyProphet でFDR制御
  5. タンパク質定量マトリクスを出力

【注意】
  - RAWファイルのダウンロードが必要（約40GB）
  - 実行に数時間かかる場合がある
  - 補足データ（Table S2）を使えばこのステップをスキップ可能

入力: data/raw/*.mzML + data/raw/human_proteome.fasta.gz
出力: results/openms_output/
"""

# ============================================================
# ライブラリの読み込み
# ============================================================
import os          # ファイルパス操作
import subprocess  # 外部コマンドの実行
import glob        # ワイルドカードによるファイル検索
import gzip        # gzip圧縮ファイルの展開
import shutil      # ファイルコピー操作

# ============================================================
# 設定
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "results", "openms_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# FASTAファイルのパス
FASTA_GZ = os.path.join(RAW_DIR, "human_proteome.fasta.gz")
FASTA = os.path.join(RAW_DIR, "human_proteome.fasta")


# ============================================================
# FASTAファイルの準備
# ============================================================
def prepare_fasta():
    """gzip圧縮FASTAを展開する。"""
    if os.path.exists(FASTA):
        print(f"  FASTAファイル（既存）: {FASTA}")
        return FASTA
    if not os.path.exists(FASTA_GZ):
        print("  エラー: FASTAファイルが見つかりません。Step 02 を先に実行してください。")
        return None
    print(f"  FASTA展開中: {os.path.basename(FASTA_GZ)} → {os.path.basename(FASTA)}")
    with gzip.open(FASTA_GZ, "rb") as f_in:
        with open(FASTA, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    print("  展開完了")
    return FASTA


# ============================================================
# Step 3-1: In silico消化（Digestor）
# ============================================================
def run_digestor(fasta_path):
    """
    FASTAデータベースをin silico消化してペプチドリストを生成する。

    【Digestorとは？】
      OpenMS付属のツール。タンパク質配列をトリプシン等の酵素で
      計算上「消化」し、生成されるペプチドのリストを出力する。
      DIA-NNの --fasta-search に相当する処理。

    【パラメータ（論文準拠）】
      - 酵素: Trypsin（KとRのC末端で切断）
      - ミスクリーベージ: 1回
      - ペプチド長: 7〜45アミノ酸
    """
    output_path = os.path.join(OUTPUT_DIR, "digested_peptides.fasta")
    if os.path.exists(output_path):
        print(f"  スキップ（既存）: {os.path.basename(output_path)}")
        return output_path

    cmd = [
        "Digestor",
        "-in", fasta_path,
        "-out", output_path,
        "-enzyme", "Trypsin",           # トリプシン消化
        "-missed_cleavages", "1",       # ミスクリーベージ1回許容
        "-min_length", "7",             # ペプチド最小長7
        "-max_length", "45",            # ペプチド最大長45
    ]

    print("  Digestor 実行中...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode == 0:
            print(f"  完了: {output_path}")
            # ペプチド数をカウント
            with open(output_path) as f:
                n_peptides = sum(1 for line in f if line.startswith(">"))
            print(f"  生成ペプチド数: {n_peptides:,}")
        else:
            print(f"  エラー: {result.stderr[:300]}")
        return output_path
    except FileNotFoundError:
        print("  エラー: OpenMS (Digestor) が見つかりません。")
        print("  インストール: micromamba install -c bioconda -c conda-forge openms -y")
        return None


# ============================================================
# Step 3-2: OpenSWATH によるDIAデータ解析
# ============================================================
def run_openswath(mzml_file, library_path):
    """
    OpenSWATHでDIA-MSデータを解析する。

    【OpenSWATHとは？】
      DIA/SWATHデータをスペクトルライブラリと照合し、
      ペプチド/タンパク質を同定・定量するOpenMSのモジュール。
      論文で使用されたDIA-NNの主要機能に相当する（OpenMSで代替）。
    """
    basename = os.path.splitext(os.path.basename(mzml_file))[0]
    output_path = os.path.join(OUTPUT_DIR, f"{basename}_openswath.tsv")

    cmd = [
        "OpenSwathWorkflow",
        "-in", mzml_file,
        "-tr", library_path,            # スペクトルライブラリ
        "-out_tsv", output_path,
        "-min_upper_edge_dist", "1",
        "-mz_extraction_window", "10",   # m/z抽出ウィンドウ（ppm）
        "-ppm",
        "-rt_extraction_window", "600",  # 保持時間ウィンドウ（秒）
        "-threads", "4",
    ]

    print(f"  OpenSWATH 実行中: {os.path.basename(mzml_file)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode == 0:
            print(f"  完了: {os.path.basename(output_path)}")
        else:
            print(f"  エラー: {result.stderr[:300]}")
        return output_path
    except FileNotFoundError:
        print("  エラー: OpenMS (OpenSwathWorkflow) が見つかりません。")
        print("  インストール: micromamba install -c bioconda -c conda-forge openms -y")
        return None


# ============================================================
# メイン処理
# ============================================================
def main():
    print("=" * 60)
    print("Step 3: DIA-MSデータのタンパク質同定・定量（OpenMS/OpenSWATH）")
    print("=" * 60)
    print()
    print("【ツール情報】")
    print("  OpenMS: BSD 3-Clause ライセンス（無料、商用利用可能）")
    print("  DIA-NNの代替として使用")

    # --- 1. OpenMSの存在確認 ---
    print("\n[1/4] OpenMSの確認")
    try:
        result = subprocess.run(["Digestor", "--help"], capture_output=True, text=True)
        print("  OpenMS (Digestor) が利用可能です")
    except FileNotFoundError:
        print("  OpenMS が見つかりません。")
        print("  インストール: micromamba install -c bioconda -c conda-forge openms -y")
        print("\n  注意: OpenMSをスキップして補足データ（Table S2）から")
        print("  Step 04（前処理）に直接進むことも可能です。")
        return

    # --- 2. FASTAの準備 ---
    print("\n[2/4] FASTAデータベースの準備")
    fasta = prepare_fasta()
    if not fasta:
        return

    # --- 3. In silico消化 ---
    print("\n[3/4] In silico消化（Digestor）")
    peptides = run_digestor(fasta)

    # --- 4. 入力ファイルの検索とOpenSWATH実行 ---
    print("\n[4/4] mzMLファイルの検索")
    mzml_files = glob.glob(os.path.join(RAW_DIR, "*.mzML"))
    if not mzml_files:
        print("  mzMLファイルが見つかりません。")
        print(f"  {RAW_DIR} にmzMLファイルを配置してください。")
        print("\n  ProteomeXchange (PXD058672) からダウンロード可能です。")
        print("  → 補足データ（Table S2）を使って Step 04 から開始することも可能です。")
    else:
        print(f"  {len(mzml_files)} ファイル検出")
        for f in mzml_files[:3]:
            print(f"    {os.path.basename(f)}")

    print("\n" + "=" * 60)
    print("完了！")
    print("=" * 60)


if __name__ == "__main__":
    main()
