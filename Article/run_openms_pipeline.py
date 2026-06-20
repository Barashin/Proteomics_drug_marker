#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenMS + AlphaPeptDeep DIA解析パイプライン
GPU環境（RTX 4070）で実行 - Toyota et al. 2025論文再現
"""

import os, glob, re, subprocess, sys, time
import pandas as pd
import numpy as np

def main():
    print("=== OpenMS + AlphaPeptDeep パイプライン開始 ===")
    print("GPU環境（RTX 4070）+ Toyota et al. 2025パラメータで実行中...")

    # --- パス設定 ---
    MZML_DIR = "data/raw/raw_mzML"
    FASTA_PATH = "data/raw/human_proteome.fasta"
    RESULTS_DIR = "results"

    # OpenMS + AlphaPeptDeep パイプライン専用ディレクトリ
    OPENMS_OUT = os.path.join(RESULTS_DIR, "openms_output")
    LIBRARY_DIR = os.path.join(OPENMS_OUT, "library")
    OSWR_DIR = os.path.join(OPENMS_OUT, "openswath")
    PYPROPHET_DIR = os.path.join(OPENMS_OUT, "pyprophet")

    # ディレクトリ作成
    for d in [OPENMS_OUT, LIBRARY_DIR, OSWR_DIR, PYPROPHET_DIR]:
        os.makedirs(d, exist_ok=True)

    # mzMLファイル一覧
    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
    print(f"mzMLファイル: {len(mzml_files)}個")

    # --- Toyota et al. 2025論文パラメータ ---
    Q_THRESHOLD = 0.01  # 1% FDR
    MISSED_CLEAVAGES = 1
    PEPTIDE_MIN_LEN = 7   # Toyota論文: 7-45
    PEPTIDE_MAX_LEN = 45
    PRECURSOR_CHARGE_MIN = 2  # Toyota論文: 2-4
    PRECURSOR_CHARGE_MAX = 4
    FRAGMENT_MIN_MZ = 200.0   # Toyota論文: 200-1800
    FRAGMENT_MAX_MZ = 1800.0
    MS1_PPM = 10.0  # Toyota論文: 10 ppm
    MS2_PPM = 10.0

    print("Toyota et al. 2025パラメータ:")
    print(f"  ペプチド長: {PEPTIDE_MIN_LEN}-{PEPTIDE_MAX_LEN}")
    print(f"  電荷: {PRECURSOR_CHARGE_MIN}-{PRECURSOR_CHARGE_MAX}")
    print(f"  フラグメント: {FRAGMENT_MIN_MZ}-{FRAGMENT_MAX_MZ}")
    print(f"  質量精度: {MS1_PPM}/{MS2_PPM} ppm")

    # === Step 1: AlphaPeptDeep スペクトルライブラリ生成 ===
    print("\n=== Step 1: AlphaPeptDeep (GPU加速) ===")

    # 既にGPU環境で生成済みのライブラリがあるかチェック
    PREDICTED_LIB = os.path.join(LIBRARY_DIR, "predicted_library_gpu.tsv")
    if os.path.exists(PREDICTED_LIB):
        lib_size = os.path.getsize(PREDICTED_LIB) / 1024**2
        print(f"既存のGPUライブラリを使用: {lib_size:.1f} MB")
        print(f"ファイル: {PREDICTED_LIB}")
    else:
        # GPU環境で生成済みのライブラリをコピー
        gpu_lib = "results/openms_output/library_gpu/predicted_library_gpu.tsv"
        if os.path.exists(gpu_lib):
            import shutil
            shutil.copy2(gpu_lib, PREDICTED_LIB)
            lib_size = os.path.getsize(PREDICTED_LIB) / 1024**2
            print(f"GPUライブラリをコピー: {lib_size:.1f} MB")
        else:
            print("エラー: GPUライブラリが見つかりません")
            return

    # ライブラリ確認
    if os.path.exists(PREDICTED_LIB):
        df_lib = pd.read_csv(PREDICTED_LIB, sep="\t", nrows=5)
        print(f"ライブラリ列: {len(df_lib.columns)}個")
        print(f"例: {list(df_lib.columns[:5])}")
    else:
        print("エラー: ライブラリファイルが見つかりません")
        return

    print("Step 1完了: AlphaPeptDeepライブラリ準備済み")

if __name__ == "__main__":
    main()