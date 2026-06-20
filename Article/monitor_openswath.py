#!/usr/bin/env python3
"""
OpenSWATH進行監視スクリプト
32サンプル全完了時に自動的にPyProphetワークフローを実行
"""

import os, glob, time, subprocess
import datetime

def check_openswath_progress():
    """OpenSWATH進行状況確認"""
    oswr_dir = "results/openms_output/openswath"
    osw_files = glob.glob(os.path.join(oswr_dir, "*.osw"))

    # ファイルサイズ確認（完了判定）
    complete_files = []
    for osw_file in osw_files:
        size = os.path.getsize(osw_file)
        if size > 590000000:  # 590MB以上なら完了とみなす
            complete_files.append(osw_file)

    return len(complete_files), len(osw_files)

def main():
    print("=== OpenSWATH監視開始 ===")
    print("32サンプル完了まで監視し、完了時にPyProphetを自動実行")

    target_samples = 32
    check_interval = 60  # 1分間隔

    while True:
        complete, total = check_openswath_progress()
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")

        print(f"[{timestamp}] 進行状況: {complete}/{target_samples} サンプル完了")

        if complete >= target_samples:
            print("🎉 全サンプル完了！PyProphetワークフロー開始")

            # PyProphetワークフロー実行
            subprocess.run(["python", "run_pyprophet_workflow.py"])
            break

        print(f"  -> {target_samples - complete}サンプル残り、{check_interval}秒後再確認")
        time.sleep(check_interval)

if __name__ == "__main__":
    main()