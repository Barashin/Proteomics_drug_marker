#!/usr/bin/env python3
"""
Step 3: mzMLファイルの中身を調べる
==================================
DIA解析を始める前に、手元にあるmzMLファイルが想定通りの中身になっているかを
確認します。具体的には以下のメタ情報を抽出し、書籍・ブログ記事の数値ソースとして使います。

  - ファイル数とサイズ
  - 装置情報（機種名、ベンダー）
  - 総スペクトル数、MS1/MS2の内訳
  - MS1/MS2 の m/z 範囲
  - 保持時間（RT）範囲
  - centroid化されているか（centroid=True/False）

【なぜこのステップが必要か？】
- rawファイルから mzML 変換する際のパラメータ（ピークピッキング等）によって、
  ダウンストリーム解析の挙動が変わるため、事前に確認しておく必要がある
- DIA解析ツール（sage等）は centroided MS2 を前提とするツールが多い
- 装置情報は論文の記述と一致しているかの確認にも使う

入力: data/raw/raw_mzML/*.mzML
出力: results/tables/mzml_inventory.csv
"""

# ============================================================
# ライブラリの読み込み
# ============================================================
import os                 # ファイルパス操作
import glob               # パターンマッチによるファイル列挙
import pandas as pd       # データフレーム操作
import pymzml             # mzMLパーサ（MIT ライセンス）

# ============================================================
# 設定（パスの定義）
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
MZML_DIR = os.path.join(PROJECT_DIR, "data", "raw", "raw_mzML")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

OUTPUT_CSV = os.path.join(TABLES_DIR, "mzml_inventory.csv")


# ============================================================
# ヘルパー関数
# ============================================================
def _is_centroided(element) -> bool:
    """spectrum 要素内の cvParam を走査し centroid spectrum (MS:1000127) を探す。

    mzMLでは各スペクトルの centroid/profile 区別は
    `<cvParam accession="MS:1000127" name="centroid spectrum" />` または
    `<cvParam accession="MS:1000128" name="profile spectrum" />` で示される。
    """
    # 名前空間を無視して accession 属性を確認する
    for child in element.iter():
        acc = child.attrib.get("accession", "")
        if acc == "MS:1000127":
            return True
        if acc == "MS:1000128":
            return False
    return False


# ============================================================
# mzML検査関数
# ============================================================
def inspect_mzml(path: str) -> dict:
    """1つのmzMLファイルを走査し、メタ情報を辞書で返す。

    pymzml.run.Reader は lazy にスペクトルをイテレートするので、
    1ファイルあたり数十秒〜1分程度で全スペクトルをパースできる。
    """
    filename = os.path.basename(path)
    filesize_mb = round(os.path.getsize(path) / (1024 * 1024), 1)
    print(f"[inspect] {filename} ({filesize_mb} MB) ...", flush=True)

    reader = pymzml.run.Reader(path)

    # 装置情報（mzML header から取得）
    instrument = ""
    try:
        # pymzml は reader.info に OBO情報を格納する
        instrument_list = reader.info.get("referenceable_param_group_list", [])
        if instrument_list:
            instrument = str(instrument_list)[:120]
    except Exception:
        pass

    total_spectra = 0
    ms1_count = 0
    ms2_count = 0
    ms1_mz_min = float("inf")
    ms1_mz_max = float("-inf")
    ms2_mz_min = float("inf")
    ms2_mz_max = float("-inf")
    rt_min = float("inf")
    rt_max = float("-inf")
    centroided_ms1 = None
    centroided_ms2 = None
    isolation_windows = set()

    for spec in reader:
        total_spectra += 1
        ms_level = spec.ms_level
        rt = spec.scan_time_in_minutes() if spec.scan_time_in_minutes() is not None else None
        if rt is not None:
            if rt < rt_min:
                rt_min = rt
            if rt > rt_max:
                rt_max = rt

        if ms_level == 1:
            ms1_count += 1
            if len(spec.peaks("raw")) > 0:
                mz_arr = spec.peaks("raw")[:, 0]
                lo, hi = float(mz_arr.min()), float(mz_arr.max())
                if lo < ms1_mz_min:
                    ms1_mz_min = lo
                if hi > ms1_mz_max:
                    ms1_mz_max = hi
            if centroided_ms1 is None:
                # MS:1000127 = "centroid spectrum", MS:1000128 = "profile spectrum"
                # cvParam 子要素の accession を走査する
                centroided_ms1 = _is_centroided(spec.element)

        elif ms_level == 2:
            ms2_count += 1
            if len(spec.peaks("raw")) > 0:
                mz_arr = spec.peaks("raw")[:, 0]
                lo, hi = float(mz_arr.min()), float(mz_arr.max())
                if lo < ms2_mz_min:
                    ms2_mz_min = lo
                if hi > ms2_mz_max:
                    ms2_mz_max = hi
            if centroided_ms2 is None:
                centroided_ms2 = _is_centroided(spec.element)

            # Isolation window 幅（DIA判定用）
            try:
                precursors = spec.selected_precursors
                if precursors:
                    # pymzml は isolation window を直接は出さないが
                    # selected_precursors[0]["mz"] で target m/z を取得可能
                    target_mz = precursors[0].get("mz")
                    if target_mz is not None:
                        isolation_windows.add(round(float(target_mz), 1))
            except Exception:
                pass

    reader.close() if hasattr(reader, "close") else None

    return {
        "file": filename,
        "size_MB": filesize_mb,
        "total_spectra": total_spectra,
        "ms1_count": ms1_count,
        "ms2_count": ms2_count,
        "ms1_mz_min": round(ms1_mz_min, 2) if ms1_mz_min != float("inf") else None,
        "ms1_mz_max": round(ms1_mz_max, 2) if ms1_mz_max != float("-inf") else None,
        "ms2_mz_min": round(ms2_mz_min, 2) if ms2_mz_min != float("inf") else None,
        "ms2_mz_max": round(ms2_mz_max, 2) if ms2_mz_max != float("-inf") else None,
        "rt_min_min": round(rt_min, 2) if rt_min != float("inf") else None,
        "rt_max_min": round(rt_max, 2) if rt_max != float("-inf") else None,
        "centroided_ms1": centroided_ms1,
        "centroided_ms2": centroided_ms2,
        "n_isolation_targets": len(isolation_windows),
    }


# ============================================================
# メイン
# ============================================================
def main():
    mzml_files = sorted(glob.glob(os.path.join(MZML_DIR, "*.mzML")))
    if not mzml_files:
        raise SystemExit(f"No mzML files found in {MZML_DIR}")
    print(f"Found {len(mzml_files)} mzML files under {MZML_DIR}")

    records = []
    for path in mzml_files:
        try:
            records.append(inspect_mzml(path))
        except Exception as e:
            print(f"  ! failed: {path}: {e}")
            records.append({"file": os.path.basename(path), "error": str(e)})

    df = pd.DataFrame(records)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nWrote {OUTPUT_CSV}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
