#!/usr/bin/env python3
"""
Step 2: データ前処理
==================
DIA解析の出力（タンパク質定量データ）を読み込み、統計解析に適した形に前処理します。
論文ではPerseus（商用利用不可）を使用していますが、本スクリプトではPythonで完全代替します。

【処理の流れ】
  1. タンパク質定量データの読み込み（Table S2 or DIA解析出力）
  2. Log2変換（強度値のスケール変換）
  3. 欠損値フィルタリング（少なくとも1群で70%以上の有効値）
  4. 欠損値補完（Perseus互換: downshift法）

【なぜ前処理が必要か？】
  - 質量分析の生データは値が非常に大きく（10^6〜10^9程度）、分布が偏っている
  - Log2変換で正規分布に近づけ、統計検定が適用しやすくなる
  - 欠損値があると多くの統計手法が使えないため、補完が必要

入力: data/raw/ (DIA解析出力 or 補足Table S2)
出力: results/preprocessed_data.csv
      results/sample_info.csv
"""

# ============================================================
# ライブラリの読み込み
# ============================================================
import os            # ファイルパス操作
import numpy as np   # 数値計算ライブラリ（log2変換、乱数生成等）
import pandas as pd  # データフレーム操作ライブラリ（表形式データの読み書き）

# ============================================================
# 設定（パスとパラメータの定義）
# ============================================================

# ディレクトリパスの設定
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- 欠損値補完パラメータ（Perseus互換） ---
# Perseus の "Replace missing values from normal distribution" と同じ設定。
# downshift: 有効値の平均から何SD下にシフトするか
#   → 2.4 = 「検出されなかったタンパク質は、有効値の平均より2.4SD低い」と仮定
# width: 補完値の分布幅（SDの倍率）
#   → 0.3 = 元のSDの30%の幅で補完値をばらつかせる
IMPUTATION_DOWNSHIFT = 2.4
IMPUTATION_WIDTH = 0.3

# --- 有効値の最低割合 ---
# 少なくとも1つの群（Normal or Tumor）で、この割合以上のサンプルで
# 有効値（非欠損値）があるタンパク質のみを残す。
# 0.70 = 70%（論文の設定に準拠）
VALID_VALUE_RATIO = 0.70


# ============================================================
# データ読み込み関数
# ============================================================
def load_protein_data(data_dir):
    """
    タンパク質定量データを読み込む関数。

    以下の優先順位で読み込みを試みる:
      1. sage 出力由来の protein matrix（results/protein_matrix_from_sage.csv）
      2. 補足データのExcelファイル（Table S2）
      3. DIA解析の出力ファイル（report.tsv）
      4. どれも無ければデモデータを自動生成

    引数:
        data_dir (str): データディレクトリのパス

    戻り値:
        pd.DataFrame: タンパク質（行）× サンプル（列）の定量マトリクス
    """
    # --- (1) sage 由来の protein matrix があれば最優先で使う ---
    # step_04_build_protein_matrix.py が results/protein_matrix_from_sage.csv を出力している
    sage_matrix_path = os.path.join(RESULTS_DIR, "protein_matrix_from_sage.csv")
    if os.path.exists(sage_matrix_path):
        print(f"  読み込み: {sage_matrix_path} (sage DIA解析由来)")
        df = pd.read_csv(sage_matrix_path, index_col=0)
        df.index.name = "Protein"
        print(f"  読み込み完了: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
        print(f"  サンプル例: {df.columns.tolist()[:4]} ...")
        return df

    # --- Table S2 (Protein expression data) を探す ---
    # Toyota et al. 2025 の Table S2 は特殊な構造:
    #   行0: カラム名（Case 1 Non-tumor, Case 1 Tumor, ...）
    #   列0: UniProt Accession Number
    #   列1: Protein Name
    #   列2: Gene Symbol
    #   列3-34: Z-score（32サンプル）← 使わない
    #   列35-66: Log2(Protein intensity)（32サンプル）← これを使う
    #   列67: Fold-change, 列68: p-value
    s2_file = [f for f in os.listdir(data_dir) if "S2" in f and f.endswith(".xlsx")]
    if s2_file:
        filepath = os.path.join(data_dir, s2_file[0])
        print(f"  読み込み: {s2_file[0]}")

        # Table S2 の構造（Excelのセル結合あり）:
        #   行0: セクションヘッダー（"Z score", "Log2(Protein intensity)"）← セル結合
        #   行1: サンプル名（"Case 1 Non-tumor tissue", "Case 1 Tumor tissue", ...）
        #   行2以降: データ
        #   列0: UniProt Accession Number
        #   列1: Protein Name
        #   列2: Gene Symbol
        #   列3-34: Z-score（32サンプル）← 使わない
        #   列35-66: Log2(Protein intensity)（32サンプル）← これを使う
        raw = pd.read_excel(filepath, sheet_name=0, header=None)

        # 行0のセクションヘッダーから Log2 列の開始位置を特定
        header_row0 = raw.iloc[0].tolist()
        log2_start = None
        for i, h in enumerate(header_row0):
            if h is not None and "Log2" in str(h):
                log2_start = i
                break

        if log2_start is not None:
            # 行1からサンプル名を取得（Log2セクションの32カラム分）
            sample_names_raw = raw.iloc[1, log2_start:log2_start + 32].tolist()

            # Gene Symbol（列2）をインデックスに使用
            gene_symbols = raw.iloc[2:, 2].tolist()  # 行2以降のデータ

            # Log2強度値を抽出（行2以降、列log2_start〜log2_start+32）
            df = raw.iloc[2:, log2_start:log2_start + 32].copy()
            df.index = gene_symbols
            df.index.name = "Protein"

            # カラム名を整理（改行を除去 → CRC01-N / CRC01-T 形式に変換）
            # 注意: Table S2 に Case 12 Non-tumor が2回出現する（Case 11 Non-tumor が欠落）
            # 重複カラム名を自動回避するため、出現回数を追跡する
            new_cols = []
            seen = {}
            for col_name in sample_names_raw:
                col_str = str(col_name).replace("\n", " ").strip()
                if "Non-tumor" in col_str or "non-tumor" in col_str.lower():
                    case_num = col_str.split("Case")[1].strip().split()[0]
                    col_id = f"CRC{case_num.zfill(2)}-N"
                elif "Tumor" in col_str:
                    case_num = col_str.split("Case")[1].strip().split()[0]
                    col_id = f"CRC{case_num.zfill(2)}-T"
                else:
                    col_id = col_str
                # 重複チェック: 既に同じ名前があれば _2 を付与
                if col_id in seen:
                    seen[col_id] += 1
                    col_id = f"{col_id}_dup{seen[col_id]}"
                else:
                    seen[col_id] = 1
                new_cols.append(col_id)
            df.columns = new_cols

            # 数値に変換（セル結合やテキストが混入している場合に備える）
            df = df.apply(pd.to_numeric, errors="coerce")

            # Gene Symbol の重複を処理
            # 同じ遺伝子名のタンパク質が複数ある場合（アイソフォーム等）、
            # 名前に _2, _3 を付与して一意にする
            dup_mask = df.index.duplicated(keep="first")
            n_dup = dup_mask.sum()
            if n_dup > 0:
                new_index = list(df.index)
                counter = {}
                for i, name in enumerate(new_index):
                    if name in counter:
                        counter[name] += 1
                        new_index[i] = f"{name}_{counter[name]}"
                    else:
                        counter[name] = 1
                df.index = new_index
                print(f"  注意: {n_dup} 個の重複Gene Symbolを一意化しました")

            print(f"  Log2強度値: 列{log2_start}〜{log2_start+31} ({len(new_cols)} サンプル)")
            print(f"  読み込み完了: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
            print(f"  サンプル例: {new_cols[:4]} ...")
            return df

        print("  注意: Log2列が見つかりませんでした。")

    # --- その他の Excel ファイルを探す ---
    xlsx_files = [f for f in os.listdir(data_dir) if f.endswith(".xlsx") and "S2" not in f]
    for xlsx in xlsx_files:
        filepath = os.path.join(data_dir, xlsx)
        try:
            xls = pd.ExcelFile(filepath)
            for sheet in xls.sheet_names:
                if "abundance" in sheet.lower() or "expression" in sheet.lower():
                    print(f"  読み込み: {xlsx} -> シート '{sheet}'")
                    return pd.read_excel(filepath, sheet_name=sheet)
        except Exception:
            continue

    # --- 実データがない場合はデモデータを生成 ---
    print("  注意: 実データが見つかりません。デモデータを生成します。")
    return generate_demo_data()


# ============================================================
# デモデータ生成関数
# ============================================================
def generate_demo_data():
    """
    論文の記述に基づいたデモデータを生成する関数。
    実データがなくてもパイプライン全体をテスト実行できるようにする。

    【生成するデータの仕様】
      - 500 タンパク質（実データは10,329だが、デモ用に縮小）
      - 16 患者 × 2条件（Normal, Tumor）= 32 サンプル
      - 最初の100タンパク質は腫瘍で発現上昇
      - 次の100タンパク質は腫瘍で発現低下
      - 約10%の欠損値を導入
    """
    # 乱数のシード固定（毎回同じデモデータを生成するため）
    np.random.seed(42)

    n_proteins = 500  # タンパク質数（デモ用に縮小）
    n_patients = 16   # 患者数（論文と同じ）

    # 各患者のステージ情報（論文Table S1に基づく）
    stages = [
        "I", "I", "I",              # Stage I: 3人
        "II", "II", "II", "II",     # Stage II: 4人
        "III", "III", "III", "III",  # Stage III: 4人
        "IV", "IV", "IV", "IV", "IV" # Stage IV: 5人
    ]

    # サンプル名の生成（CRC01-N = 患者1の正常組織, CRC01-T = 患者1の腫瘍組織）
    samples_n = [f"CRC{i+1:02d}-N" for i in range(n_patients)]  # Normal
    samples_t = [f"CRC{i+1:02d}-T" for i in range(n_patients)]  # Tumor
    all_samples = samples_n + samples_t  # 全32サンプル

    # タンパク質名の生成
    proteins = [f"Protein_{i+1:04d}" for i in range(n_proteins)]

    # 基本発現量（log2スケール、平均20、標準偏差3の正規分布）
    # 実データのlog2変換後の値は15〜25程度
    base_expr = np.random.normal(20, 3, n_proteins)

    # 各サンプルの発現データを生成
    data = {}
    for i, sample in enumerate(all_samples):
        is_tumor = i >= n_patients  # 後半16サンプルが腫瘍

        # ノイズ（生物学的・技術的変動をシミュレート）
        noise = np.random.normal(0, 0.5, n_proteins)

        if is_tumor:
            # 腫瘍組織では一部のタンパク質が変動
            change = np.zeros(n_proteins)
            change[:100] = np.random.uniform(1, 3, 100)     # 最初の100個: 上昇
            change[100:200] = np.random.uniform(-3, -1, 100) # 次の100個: 低下
            expr = base_expr + noise + change
        else:
            expr = base_expr + noise

        # 約10%の値を欠損（NaN）にする
        # 質量分析では低発現タンパク質が検出されないことがある
        mask = np.random.random(n_proteins) < 0.10
        expr[mask] = np.nan
        data[sample] = expr

    # DataFrameに変換（行=タンパク質、列=サンプル）
    df = pd.DataFrame(data, index=proteins)
    df.index.name = "Protein"

    # 臨床情報も保存（ステージ別解析で使用）
    clinical = pd.DataFrame({
        "Sample_N": samples_n,
        "Sample_T": samples_t,
        "Stage": stages
    })
    clinical.to_csv(os.path.join(RESULTS_DIR, "clinical_info.csv"), index=False)

    print(f"  デモデータ生成: {n_proteins} タンパク質 × {len(all_samples)} サンプル")
    return df


# ============================================================
# Log2変換関数
# ============================================================
def log2_transform(df):
    """
    生の強度値をLog2変換する関数。
    既にLog2スケールの場合は自動的にスキップする。

    【なぜLog2変換するか？】
      - 質量分析の強度値は桁が非常に大きい（10^6〜10^9）
      - そのままでは統計検定の仮定（正規分布）を満たさない
      - Log2変換することで正規分布に近づき、統計手法が適用しやすくなる
      - Log2スケールでは差が「倍率」に対応（差1.0 = 2倍の変化）

    引数:
        df (pd.DataFrame): タンパク質定量データ

    戻り値:
        pd.DataFrame: Log2変換後のデータ
    """
    # 値の範囲で変換済みかどうかを判定
    # Log2変換済み: 中央値が15〜25程度
    # 未変換（生データ）: 中央値が10^3以上
    median_val = df.median().median()

    if median_val > 100:
        print(f"  Log2変換実行（中央値={median_val:.1f} → 生データと判定）")
        # 0の値はLog2(0)=-∞になるため、先にNaNに置換してから変換
        df = np.log2(df.replace(0, np.nan))
    else:
        print(f"  Log2変換スキップ（中央値={median_val:.1f} → 既にlog2スケールと判定）")

    return df


# ============================================================
# 欠損値フィルタリング関数
# ============================================================
def filter_valid_values(df, sample_groups, ratio=VALID_VALUE_RATIO):
    """
    有効値が少なすぎるタンパク質を除去する関数。
    「少なくとも1つの群で、指定割合以上のサンプルに有効値がある」タンパク質のみ残す。

    【なぜフィルタリングするか？】
      - ほとんどのサンプルで検出されないタンパク質は、統計的に信頼できない
      - 欠損値が多すぎると、補完の精度も低下する
      - 論文では「少なくとも1群で70%以上」という基準を使用

    引数:
        df (pd.DataFrame): タンパク質定量データ
        sample_groups (dict): 群名→サンプル名リストの辞書
        ratio (float): 有効値の最低割合（デフォルト: 0.70）

    戻り値:
        pd.DataFrame: フィルタリング後のデータ
    """
    n_before = len(df)
    # 各タンパク質について、いずれかの群で基準を満たすかチェック
    keep_mask = pd.Series(False, index=df.index)

    for group_name, samples in sample_groups.items():
        group_data = df[samples]
        # 各タンパク質の有効値割合を計算
        # notna() → True/Falseの行列、sum(axis=1) → 各行のTrue数
        valid_ratio = group_data.notna().sum(axis=1) / len(samples)
        # いずれかの群で基準を満たせばTrue
        keep_mask |= (valid_ratio >= ratio)

    df_filtered = df[keep_mask]
    n_after = len(df_filtered)
    print(f"  有効値フィルタ: {n_before} → {n_after} タンパク質（{n_before - n_after} 除去）")
    return df_filtered


# ============================================================
# 欠損値補完関数（Perseus互換 downshift法）
# ============================================================
def impute_missing_values(df, downshift=IMPUTATION_DOWNSHIFT, width=IMPUTATION_WIDTH):
    """
    Perseus互換の欠損値補完を行う関数。

    【補完の考え方】
      質量分析で「検出されなかった」タンパク質は、存在しないのではなく
      「発現量が低すぎて検出限界以下だった」と考える。
      そこで、検出された値の分布から推測して低い値を割り当てる。

    【具体的な手順（各サンプルごとに）】
      1. 有効値（検出された値）の平均と標準偏差を計算
      2. 補完値の中心 = 平均 - downshift × 標準偏差
         → 有効値の分布より downshift × SD 分だけ低い位置
      3. 補完値の幅 = width × 標準偏差
         → 元の分布より狭い（すべて低い値に集中させるため）
      4. この分布から欠損値の数だけランダムにサンプリング

    引数:
        df (pd.DataFrame): タンパク質定量データ（欠損値含む）
        downshift (float): 平均からのシフト量（SD単位）
        width (float): 補完値の分布幅（SDの倍率）

    戻り値:
        pd.DataFrame: 欠損値が補完されたデータ
    """
    df_imputed = df.copy()  # 元データを変更しないようにコピー
    n_imputed = 0  # 補完した値の総数カウンタ

    # 各サンプル（列）ごとに補完
    for col in df.columns:
        # そのサンプルの有効値（NaNでない値）を取得
        valid = df[col].dropna()
        if len(valid) == 0:
            continue  # 全欠損の場合はスキップ

        # 有効値の統計量
        col_mean = valid.mean()  # 平均
        col_std = valid.std()    # 標準偏差

        # 補完値の分布パラメータを計算
        # 例: 平均=20, SD=2, downshift=2.4, width=0.3 の場合
        #   imp_mean = 20 - 2.4 × 2 = 15.2（有効値の平均より4.8低い）
        #   imp_std = 0.3 × 2 = 0.6（元のSDの30%の幅）
        imp_mean = col_mean - downshift * col_std
        imp_std = width * col_std

        # 欠損値の位置を特定（iloc を使ってインデックス重複問題を回避）
        col_idx = df.columns.get_loc(col)
        col_values = df_imputed.iloc[:, col_idx]
        missing_mask = col_values.isna()
        n_missing = int(missing_mask.sum())

        if n_missing > 0:
            # 正規分布からランダムにサンプリングして補完
            imp_values = np.random.normal(imp_mean, imp_std, n_missing)
            df_imputed.iloc[missing_mask.values, col_idx] = imp_values
            n_imputed += n_missing

    print(f"  欠損値補完: {n_imputed} 値を補完（downshift={downshift}, width={width}）")
    return df_imputed


# ============================================================
# メイン処理
# ============================================================
def main():
    print("=" * 60)
    print("Step 2: データ前処理")
    print("=" * 60)

    # ----------------------------------------------------------
    # [1/4] データ読み込み
    # ----------------------------------------------------------
    print("\n[1/4] タンパク質定量データの読み込み")
    df = load_protein_data(RAW_DIR)

    print(f"\n  データサイズ: {df.shape[0]} 行 × {df.shape[1]} 列")

    # タンパク質名をインデックス（行名）に設定
    if df.index.name != "Protein" and "Protein" in df.columns:
        df = df.set_index("Protein")
    elif df.index.dtype == "int64":
        df = df.set_index(df.columns[0])

    # 数値列のみ抽出（文字列列を除外）
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    df = df[numeric_cols]

    # ----------------------------------------------------------
    # [2/4] サンプル群の定義（Normal群 と Tumor群）
    # ----------------------------------------------------------
    print("\n[2/4] サンプル群の定義")

    # サンプル名から自動的にNormal/Tumorを判別
    # Toyota et al. では "-N" = Normal, "-T" = Tumor
    normal_samples = [c for c in df.columns if "-N" in c or "Normal" in c or "non" in c.lower()]
    tumor_samples = [c for c in df.columns if "-T" in c or "Tumor" in c
                     and c not in normal_samples]

    # 自動検出できなかった場合のフォールバック
    if not normal_samples or not tumor_samples:
        mid = len(df.columns) // 2
        normal_samples = list(df.columns[:mid])
        tumor_samples = list(df.columns[mid:])

    print(f"  Normal群: {len(normal_samples)} サンプル")
    print(f"  Tumor群:  {len(tumor_samples)} サンプル")

    sample_groups = {"Normal": normal_samples, "Tumor": tumor_samples}

    # ----------------------------------------------------------
    # [3/4] Log2変換
    # ----------------------------------------------------------
    print("\n[3/4] Log2変換")
    df = log2_transform(df)

    # ----------------------------------------------------------
    # [4/4] フィルタリングと欠損値補完
    # ----------------------------------------------------------
    print("\n[4/4] 有効値フィルタリングと欠損値補完")
    df = filter_valid_values(df, sample_groups)
    df = impute_missing_values(df)

    # ----------------------------------------------------------
    # 結果の保存
    # ----------------------------------------------------------
    # 前処理済みデータを CSV で保存
    output_path = os.path.join(RESULTS_DIR, "preprocessed_data.csv")
    df.to_csv(output_path)
    print(f"\n保存: {output_path}")

    # サンプル情報を CSV で保存（後のステップで使用）
    sample_info = pd.DataFrame({
        "Sample": normal_samples + tumor_samples,
        "Condition": ["Normal"] * len(normal_samples) + ["Tumor"] * len(tumor_samples)
    })
    sample_info_path = os.path.join(RESULTS_DIR, "sample_info.csv")
    sample_info.to_csv(sample_info_path, index=False)
    print(f"保存: {sample_info_path}")

    print("\n" + "=" * 60)
    print(f"前処理完了: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
    print("=" * 60)


# ============================================================
# スクリプトの実行エントリーポイント
# ============================================================
if __name__ == "__main__":
    main()
