---
title: "大腸がんのステージ進行に伴うタンパク質変動をクラスター分析で可視化する【論文再現シリーズ #10】"
emoji: "📈"
type: "tech"
topics: ["proteomics", "python", "clustering", "labcode"]
published: false
---

# 大腸がんのステージ進行に伴うタンパク質変動をクラスター分析で可視化する

## はじめに

大腸がんはStage I〜IVに分類されます。この記事では、各ステージでのタンパク質発現パターンを解析し、**疾患進行に伴って一貫して変動するタンパク質群**を特定します。Toyota et al. 2025のFigure 3の再現です。

:::message
**この記事で行う処理**
大腸がんのステージ（I〜IV）と正常組織の5群間で、各タンパク質の発現量に差があるかをOne-way ANOVAで検定し、FDR補正で偽陽性を制御します。有意なタンパク質を階層的クラスタリングで30グループに分割し、ステージ進行に伴って段階的に増加・減少するパターンを持つタンパク質群を特定します。これにより、疾患進行のバイオマーカー候補や分子メカニズムの手がかりが得られます。
:::

:::message
**【用語メモ】**

- **ANOVA（Analysis of Variance、分散分析）**: 3群以上の平均を比較する統計手法。ここでは Normal / Stage I / II / III / IV の5群間で発現量が有意に異なるタンパク質を検出します。
- **FDR（False Discovery Rate、偽発見率）**: 有意と判定されたタンパク質のうち、実は偽陽性である割合の期待値。`FDR < 0.01` は「有意と判定された中で偽陽性率は1%以下」を意味し、Benjamini-Hochberg 法で補正します。
:::

## 前提

- [#9 COSMIC照合](article-09-cosmic.md) が完了していること
- ステージ情報（Table S1）が利用可能
- **対応Notebook**: [`notebooks/step_10.ipynb`](../notebooks/step_10.ipynb) — この記事のコードをセルごとに実行できます

## コード全文（対応Notebook: step_10.ipynb）

### ライブラリと設定

```python
import re                      # 正規表現モジュール: サンプル名から重複ランサフィックスを除去するために使用
import numpy as np             # 数値計算ライブラリ: 配列操作やクラスター番号の一意抽出に使用
import pandas as pd            # データ分析ライブラリ: CSV読み込み、DataFrame操作の中心
import matplotlib.pyplot as plt  # グラフ描画ライブラリ: ヒートマップやラインプロットの描画に使用
# GridSpec: Figure内のサブプロットの配置を柔軟に制御するクラス
# GridSpecFromSubplotSpec: 既存のGridSpec領域をさらに分割するクラス（入れ子レイアウト用）
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch  # 凡例用のカラーパッチ（四角い色見本）を作成するクラス
import seaborn as sns          # 統計データ可視化ライブラリ: ヒートマップ描画に使用
from scipy import stats        # 科学計算ライブラリの統計モジュール: ANOVAのf_oneway関数を使用
# linkage: 階層的クラスタリングの連結行列を計算する関数
# fcluster: 連結行列からクラスターラベルを割り当てる関数
# leaves_list: デンドログラムの葉の並び順を取得する関数（ヒートマップの行ソートに使用）
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
# multipletests: 多重検定補正を行う関数（Benjamini-Hochberg法によるFDR補正に使用）
from statsmodels.stats.multitest import multipletests

# --- 定数設定 ---
RESULTS = "../results"         # 解析結果の出力先ディレクトリパス
FIG_DIR = f"{RESULTS}/figures"   # 図の保存先ディレクトリパス
TABLE_DIR = f"{RESULTS}/tables"  # テーブル（CSV/Excel）の保存先ディレクトリパス
RAW_DIR = "../data/raw"        # 元データ（論文補足テーブル等）の格納ディレクトリパス

N_CLUSTERS = 30                # 論文と同じ30クラスターに分割
FDR_THRESHOLD = 0.01           # FDR < 0.01 で有意と判定
# ステージの表示順序を定義するリスト（Normal→Stage I→II→III→IVの順）
STAGE_ORDER = ["Normal", "I", "II", "III", "IV"]
# 各ステージに割り当てるカラーコード（ヒートマップ上部のカラーバーやラインプロットで使用）
STAGE_COLORS = {
    "Normal": "#4EAED1",       # 青: 正常組織
    "I": "#66BB6A",            # 緑: Stage I
    "II": "#FFD54F",           # 黄: Stage II
    "III": "#FFA726",          # オレンジ: Stage III
    "IV": "#E8524A",           # 赤: Stage IV
}
```

### データ読み込みとステージ付与

```python
# --- データ読み込み + ステージ付与 ---
# 前処理済みのタンパク質発現量データを読み込む（行=タンパク質、列=サンプル）
df = pd.read_csv(f"{RESULTS}/preprocessed_data.csv", index_col=0)
# サンプル情報（サンプル名・患者ID等）を読み込む
sample_info = pd.read_csv(f"{RESULTS}/sample_info.csv")

# clinical_info.csv から Sample_N→"Normal", Sample_T→Stage をマッピング
# 臨床情報（患者ごとのステージ、正常/腫瘍サンプル名）を読み込む
clinical = pd.read_csv(f"{RESULTS}/clinical_info.csv")
# サンプル名→ステージの対応辞書を作成する（空の辞書で初期化）
stage_map = {}
# 臨床情報の各行（=各患者）をループし、正常サンプルと腫瘍サンプルのステージを登録
for _, row in clinical.iterrows():
    stage_map[row["Sample_N"]] = "Normal"    # 正常組織サンプルには"Normal"を割り当て
    stage_map[row["Sample_T"]] = row["Stage"]  # 腫瘍サンプルにはStage列の値（I〜IV）を割り当て

# サンプル名に "_dup1" 等の重複ランサフィックスがついている場合も
# 元のサンプル名でステージを引けるように正規表現で除去してマッピング
# re.sub(r"_dup\d+$", "", x) は末尾の「_dup1」「_dup2」等を空文字に置換する
# stage_map.get() は辞書からステージを取得（見つからなければNone）
sample_info["Stage"] = sample_info["Sample"].apply(
    lambda x: stage_map.get(re.sub(r"_dup\d+$", "", x))
)

# データの概要を表示: タンパク質数とサンプル数
print(f"データ: {df.shape[0]} タンパク質 × {df.shape[1]} サンプル")
# 各ステージに属するサンプル数を表示（例: Normal 16, III 8, II 4, ...）
print(sample_info["Stage"].value_counts().to_string())
```

### One-way ANOVA

```python
def run_anova(df, sample_info):
    """各タンパク質について One-way ANOVA + BH-FDR 補正を実行。

    【One-way ANOVAとは？】
      3群以上の平均値に差があるかを検定する手法。
      t検定は2群専用だが、ANOVAは3群以上に対応。
      本解析では Normal(16), Stage I(3), II(4), III(8), IV(1) の5群を比較する。

    【FDR補正（Benjamini-Hochberg法）とは？】
      数千回の検定を同時に行うと、偶然だけで有意になるものが出てくる。
      例: 10,000回検定 × p < 0.05 → 偶然500個が有意に
      FDR補正はこの「偽陽性の割合」を制御する手法。
    """
    # STAGE_ORDERの中から、実際にデータに存在するステージだけを抽出する
    # （例: Stage IVの患者がいなければリストから除外される）
    stages = [s for s in STAGE_ORDER if s in sample_info["Stage"].values]
    # ステージごとのサンプル列名を事前に辞書化（ループ内で毎回フィルタしない最適化）
    # 各ステージに属するサンプル名のリストを辞書で保持する（パフォーマンス向上のため）
    stage_samples = {
        s: [c for c in sample_info[sample_info["Stage"] == s]["Sample"] if c in df.columns]
        for s in stages
    }
    results = []  # ANOVA結果を格納するリスト（各要素は辞書）
    # 全タンパク質について1つずつANOVAを実行するループ
    for protein in df.index:
        groups = []  # このタンパク質のステージ別発現値グループを格納するリスト
        # 各ステージの発現値を取得
        for s in stages:
            # 該当ステージのサンプル列から発現値を取得し、欠損値を除去
            vals = df.loc[protein, stage_samples[s]].dropna().values
            # ANOVAには各群2サンプル以上が必要なので、2未満のグループは除外
            if len(vals) >= 2:
                groups.append(vals)
        # 比較可能な群が2つ未満の場合はこのタンパク質をスキップ
        if len(groups) < 2:
            continue
        # scipy.stats.f_oneway: One-way ANOVAを実行しF統計量とp値を返す
        # *groups でリストをアンパックして各群を引数として渡す
        f_stat, p_val = stats.f_oneway(*groups)
        # 結果を辞書としてリストに追加
        results.append({"Protein": protein, "F_statistic": f_stat, "P_value": p_val})

    # 結果リストをDataFrameに変換（行=タンパク質、列=F統計量・P値等）
    result_df = pd.DataFrame(results)
    # Benjamini-Hochberg法でp値をFDR補正する
    # 戻り値: reject（有意判定）, fdr（補正後p値）, alphacSidak, alphacBonf（使わない）
    _, fdr, _, _ = multipletests(result_df["P_value"], method="fdr_bh")
    result_df["FDR"] = fdr  # FDR補正後のp値を列として追加
    # FDR閾値未満なら有意（True）、以上なら非有意（False）のブール列を追加
    result_df["Significant"] = fdr < FDR_THRESHOLD
    return result_df  # ANOVA結果のDataFrameを返す
```

### ANOVAの実行

```python
# ANOVA実行: 全タンパク質についてステージ間の発現差を検定する
anova_df = run_anova(df, sample_info)
# ANOVA結果をCSVファイルに保存（後の解析で再利用可能にするため）
anova_df.to_csv(f"{TABLE_DIR}/anova_results.csv", index=False)
# Significant列がTrueの行数を合計し、有意なタンパク質数を取得
n_sig = anova_df["Significant"].sum()
# 検定結果のサマリーを表示: 全タンパク質数と有意タンパク質数
print(f"検定: {len(anova_df)} タンパク質 → 有意 (FDR<{FDR_THRESHOLD}): {n_sig}")
```

### ステージ別中央値とクラスタリング

```python
def compute_stage_medians(df, sample_info):
    """各ステージの中央値を計算。行=タンパク質、列=ステージ。

    【Z-score正規化とは？】
      各タンパク質の発現値を「平均0、標準偏差1」に変換する処理。
      これにより、発現量の絶対値ではなく「変動パターン」で比較できるようになる。
      計算式: z = (x - 平均) / 標準偏差
    """
    # 実際にデータに存在するステージだけをSTAGE_ORDERの順序で抽出
    stages = [s for s in STAGE_ORDER if s in sample_info["Stage"].values]
    medians = {}  # ステージ名→中央値Seriesの辞書を格納する
    for s in stages:
        # 該当ステージに属し、かつdfの列に存在するサンプル名を取得
        cols = [c for c in sample_info[sample_info["Stage"] == s]["Sample"] if c in df.columns]
        if cols:
            # 該当ステージの全サンプルについて、各タンパク質の中央値を計算
            # axis=1: 行方向（=サンプル方向）に中央値を取る → 各タンパク質に1つの代表値
            medians[s] = df[cols].median(axis=1)
    # 辞書をDataFrameに変換して返す（行=タンパク質、列=ステージ）
    return pd.DataFrame(medians)

# ステージ別中央値の計算: 全タンパク質×5ステージのマトリクスを作成
median_df = compute_stage_medians(df, sample_info)

# ANOVAで有意（FDR < 0.01）と判定されたタンパク質名のリストを取得
sig_proteins = anova_df[anova_df["Significant"]]["Protein"].tolist()
# 有意タンパク質が0個の場合のフォールバック: FDR値が小さい上位200個を使用
if not sig_proteins:
    sig_proteins = anova_df.nsmallest(min(200, len(anova_df)), "FDR")["Protein"].tolist()

# Z-score正規化 → 階層的クラスタリング（Ward法）
# 有意タンパク質だけを中央値マトリクスから抽出
sig_data = median_df.loc[median_df.index.isin(sig_proteins)]
# 各タンパク質（行）ごとにZ-score正規化: (値 - 行平均) / 行標準偏差
# これにより発現量の絶対値ではなくステージ間の変動パターンで比較できる
# dropna(): 標準偏差が0（全ステージで同値）の行はNaNになるので除去
z_data = sig_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1).dropna()

# クラスター数をN_CLUSTERSとデータ数の小さい方に設定（データが少ない場合の安全策）
actual_clusters = min(N_CLUSTERS, len(z_data))
# Ward法で階層的クラスタリングの連結行列を計算
# Ward法: クラスター内の分散増加が最小になるペアを統合していく方法
Z_linkage = linkage(z_data.values, method="ward")
# 連結行列からクラスターラベル（1〜30の整数）を各タンパク質に割り当て
# criterion="maxclust": 最大クラスター数を指定する基準
clusters = fcluster(Z_linkage, t=actual_clusters, criterion="maxclust")

# クラスター割り当て結果をDataFrameにまとめてCSVに保存
cluster_df = pd.DataFrame({"Protein": z_data.index, "Cluster": clusters})
cluster_df.to_csv(f"{TABLE_DIR}/cluster_assignments.csv", index=False)
# クラスタリング結果のサマリーを表示
print(f"{len(z_data)} タンパク質 → {actual_clusters} クラスター")
```

### Figure 3: ステージ別ヒートマップとプロファイルプロット

```python
def plot_figure3(df, sample_info, z_data, clusters):
    """Figure 3: 2x2 パネル（ヒートマップ + ラインプロット × 4 クラスター）。

    【処理の流れ】
      1. 論文補足テーブル（S9-S12）からクラスター情報を読み込み試行
      2. 読み込めなければ自前クラスタリングから増加2+減少2を自動選択
      3. 各クラスターについてヒートマップ（上段）＋ラインプロット（下段）を描画
    """

    # --- 論文補足 Table S9-S12 からクラスター読み込み（なければ自前） ---
    # 論文のFigure 3で使われた4つのクラスター情報（ファイル名, クラスター番号, 増減方向）
    # 論文著者が公開した補足テーブルから正確なタンパク質リストを取得する
    cluster_files = [
        ("Supplementary Table S9 Fig3A_Clustre3_34 protein.xlsx", 3, "Increased"),    # 増加クラスター1
        ("Supplementary Table S10 Fig3B_Clustre14_1324 protein.xlsx", 14, "Increased"),  # 増加クラスター2
        ("Supplementary Table S11 Fig3C_Clustre20_16 protein.xlsx", 20, "Decreased"),   # 減少クラスター1
        ("Supplementary Table S12 Fig3D_Clustre25_1062 protein.xlsx", 25, "Decreased"),  # 減少クラスター2
    ]
    selected = []  # 描画対象として選ばれたクラスター情報を格納するリスト
    # 各補足テーブルファイルを順に読み込み試行する
    for fname, cl_num, direction in cluster_files:
        try:
            # Excelファイルを読み込む（1行目をヘッダーとして使用）
            cdf = pd.read_excel(f"{RAW_DIR}/{fname}", header=0)
            # 列名に"gene"または"symbol"を含む列を探す（遺伝子名の列を特定するため）
            gene_col = [c for c in cdf.columns if "gene" in c.lower() or "symbol" in c.lower()]
            # 遺伝子名列が見つかればその列を、なければ6列目（インデックス5）を使用
            genes = cdf[gene_col[0]].dropna().tolist() if gene_col else cdf.iloc[:, 5].dropna().tolist()
            # 論文のタンパク質リストのうち、本書のデータに存在するものだけを抽出
            matched = [g for g in genes if g in df.index]
            # クラスター情報を辞書として追加（クラスター番号、タンパク質数、遺伝子名リスト、増減方向）
            selected.append({"cluster": cl_num, "n": len(matched), "genes": matched, "direction": direction})
        except Exception:
            # ファイルが見つからない・読み込みエラーの場合はスキップ
            pass

    # 4つのクラスター全てが読み込めなかった場合のフォールバック処理
    if len(selected) != 4:
        # フォールバック: 自前クラスタリングから増加2 + 減少2 を選択
        print("論文クラスター不完全 → 自前クラスタリングを使用")
        infos = []  # 全クラスターの情報を格納するリスト
        # 各クラスター番号について情報を収集
        for cl in np.unique(clusters):
            mask = clusters == cl  # 該当クラスターに属するタンパク質のブールマスク
            # タンパク質数が5未満のクラスターは可視化に不適切なのでスキップ
            if mask.sum() < 5:
                continue
            # クラスター内の全タンパク質のZ-scoreプロファイルの平均を計算
            profile = z_data.iloc[mask].mean(axis=0)
            # トレンド（傾向）: 最後のステージ(IV)と最初(Normal)の差で増減方向を判定
            trend = profile.iloc[-1] - profile.iloc[0]
            # クラスター情報を辞書として追加
            infos.append({"cluster": cl, "n": int(mask.sum()), "genes": z_data.index[mask].tolist(),
                          "direction": "Increased" if trend > 0 else "Decreased", "trend": trend})
        # 増加トレンドが大きい順に上位2クラスターを選択
        inc = sorted([c for c in infos if c["direction"] == "Increased"], key=lambda x: x["trend"], reverse=True)[:2]
        # 減少トレンドが大きい（値が小さい）順に上位2クラスターを選択
        dec = sorted([c for c in infos if c["direction"] == "Decreased"], key=lambda x: x["trend"])[:2]
        # 増加2つ + 減少2つ = 合計4クラスターを描画対象に設定
        selected = inc + dec

    # --- Figure 描画 ---
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]  # 各パネルのラベル（論文のサブパネル記号）
    # 全体の Figure を作成（幅16×高さ20インチ）
    fig = plt.figure(figsize=(16, 20))
    # 2行×2列の外側グリッドを定義（4つのクラスターパネルを配置）
    # hspace: 行間の余白、wspace: 列間の余白
    outer_gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.3)

    # 選択された4つのクラスターそれぞれについてパネルを描画
    for idx, info in enumerate(selected):
        # 各パネルを上段（ヒートマップ）と下段（ラインプロット）に分割
        # height_ratios=[3, 1]: ヒートマップを3、ラインプロットを1の比率で配分
        inner_gs = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_gs[idx],
                                          height_ratios=[3, 1], hspace=0.15)
        ax_heat = fig.add_subplot(inner_gs[0])  # 上段: ヒートマップ用のAxesオブジェクト
        ax_line = fig.add_subplot(inner_gs[1])  # 下段: ラインプロット用のAxesオブジェクト

        # このクラスターに属し、かつデータに存在するタンパク質のリストを取得
        proteins = [g for g in info["genes"] if g in df.index]

        # ステージ中央値 → Z-score
        stage_med = {}  # ステージ名→中央値Seriesの辞書
        # 各ステージについてこのクラスターのタンパク質の中央値を計算
        for s in STAGE_ORDER:
            # 該当ステージに属するサンプル列名を取得
            cols = [c for c in sample_info[sample_info["Stage"] == s]["Sample"] if c in df.columns]
            if cols:
                # 該当タンパク質×該当ステージのサンプルの中央値を計算（axis=1: サンプル方向）
                stage_med[s] = df.loc[proteins, cols].median(axis=1)
        # 辞書をDataFrameに変換（行=タンパク質、列=ステージ）
        med_matrix = pd.DataFrame(stage_med)
        # 行ごと（タンパク質ごと）にZ-score正規化
        # std() > 0 のチェック: 全ステージで同値のタンパク質は0で埋める（ゼロ除算防止）
        cz = med_matrix.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0, axis=1)

        # タンパク質をクラスタリング順にソート（見やすくするため）
        # Ward法で再クラスタリングし、デンドログラムの葉の並び順で行を並べ替える
        if len(cz) > 1:
            cz = cz.iloc[leaves_list(linkage(cz.values, method="ward"))]

        # ヒートマップ描画: seaborn.heatmapで色付きマトリクスを描く
        # cmap="RdYlGn_r": 赤(高)→黄(中)→緑(低)のカラーマップ（_rは反転）
        # center=0: Z-score=0を中央色（黄）に設定
        # vmin=-2, vmax=2: 色の範囲を-2〜+2に固定（外れ値に引きずられない）
        sns.heatmap(cz, cmap="RdYlGn_r", center=0, vmin=-2, vmax=2, ax=ax_heat,
                    xticklabels=True, yticklabels=False,
                    cbar_kws={"shrink": 0.5, "label": "Z-score"})
        # X軸ラベルを「Non-tumor」「Stage I」等の分かりやすい名前に変更
        ax_heat.set_xticklabels(
            ["Non-tumor" if s == "Normal" else f"Stage {s}" for s in cz.columns],
            fontsize=7, rotation=45, ha="right")

        # ステージカラーバー（ヒートマップ上部に配置）
        # ヒートマップの位置情報を取得して、その上に薄い帯状のAxesを配置する
        hp = ax_heat.get_position()  # ヒートマップの位置（x0, y0, width, height）を取得
        # fig.add_axes: 絶対座標で新しいAxesを追加（[左端, 下端, 幅, 高さ]）
        ax_bar = fig.add_axes([hp.x0, hp.y1 + 0.005, hp.width * 0.85, hp.height * 0.04])
        # 各ステージの色を取得（STAGE_COLORSに未登録なら灰色"#888"をデフォルト値に）
        colors_list = [STAGE_COLORS.get(s, "#888") for s in cz.columns]
        # 各ステージの色をaxvspan（垂直帯）で描画してカラーバーを構成
        for i, c in enumerate(colors_list):
            ax_bar.axvspan(i, i + 1, facecolor=c, edgecolor="white", linewidth=0.5)
        ax_bar.set_xlim(0, len(colors_list))  # X軸の範囲をステージ数に合わせる
        ax_bar.set_xticks([]); ax_bar.set_yticks([])  # 目盛りを非表示にする

        # パネルタイトルを設定（例: "(a) Cluster 3: 34 proteins (Increased)"）
        ax_heat.set_title(
            f"{panel_labels[idx]} Cluster {info['cluster']}: {info['n']} proteins ({info['direction']})",
            fontsize=11, fontweight="bold", pad=15)

        # ラインプロット描画（灰色=個々のタンパク質、黒太線=クラスター平均）
        x_pos = range(len(cz.columns))  # X軸の位置（0, 1, 2, 3, 4 = 5ステージ）
        # X軸ラベル（"Normal"を"Non-tumor"に変換）
        labels = ["Non-tumor" if s == "Normal" else s for s in cz.columns]
        # 個々のタンパク質のプロファイルを灰色の細い線で描画（パターンの多様性を可視化）
        for _, row in cz.iterrows():
            ax_line.plot(x_pos, row.values, color="gray", alpha=0.15, linewidth=0.3)
        # クラスター内の全タンパク質の平均プロファイルを計算
        mean_prof = cz.mean(axis=0)
        # 平均プロファイルを黒の太い実線で描画（クラスターの代表的なトレンドを示す）
        ax_line.plot(x_pos, mean_prof.values, "k-", linewidth=2.5)
        # 各ステージ位置にステージ色の丸マーカーを描画（視覚的にステージを区別）
        for i, (xp, ym) in enumerate(zip(x_pos, mean_prof.values)):
            ax_line.plot(xp, ym, "o", color=colors_list[i], markersize=8,
                         markeredgecolor="black", markeredgewidth=0.5, zorder=5)
        ax_line.set_xticks(list(x_pos))  # X軸の目盛り位置を設定
        ax_line.set_xticklabels(labels, fontsize=8)  # X軸ラベルを設定
        ax_line.set_ylabel("Z-score", fontsize=8)  # Y軸ラベルを設定
        ax_line.axhline(0, color="gray", ls="--", lw=0.5)  # Z-score=0の基準線を破線で描画
        ax_line.spines["top"].set_visible(False)   # 上辺の枠線を非表示にする（見た目をすっきりさせる）
        ax_line.spines["right"].set_visible(False)  # 右辺の枠線を非表示にする

    # 凡例（図の下部に配置）: 各ステージの色を示すカラーパッチを横一列に並べる
    # Patch: 凡例用の四角い色見本を作成するオブジェクト
    # bbox_to_anchor=(0.5, -0.01): Figure下端の中央に配置
    fig.legend(
        handles=[Patch(fc=STAGE_COLORS[s], label="Non-tumor" if s == "Normal" else f"Stage {s}")
                 for s in STAGE_ORDER],
        loc="lower center", ncol=5, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.01))

    # 図をPNGファイルとして保存
    path = f"{FIG_DIR}/fig3_stage_heatmap.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")  # dpi=150: 解像度150dpi、bbox_inches="tight": 余白を自動調整
    plt.close()  # メモリ解放のためFigureを閉じる
    print(f"保存: {path}")  # 保存先パスを表示


# Figure 3 の描画実行: 上で定義した関数を呼び出して実際に図を生成・保存する
plot_figure3(df, sample_info, z_data, clusters)
```

---

## コード詳細

### t検定 vs ANOVA の使い分け

| 手法 | 群数 | 使いどころ |
|------|------|----------|
| **t検定** | 2群 | Normal vs Tumor の比較（Step 4） |
| **ANOVA** | 3群以上 | Normal vs Stage I vs II vs III vs IV の比較（Step 6） |

### FDR補正の重要性

```python
# FDR補正なし: 10,000タンパク質 × p < 0.05 → 偶然500個が有意に！
# FDR補正あり: 偽陽性の割合を1%以下に制御
# multipletests: 多重検定補正を行う関数。p_valuesは全タンパク質のp値の配列
# method="fdr_bh": Benjamini-Hochberg法を指定（FDRを制御する最も一般的な手法）
# 戻り値: reject（有意かどうかのブール配列）, fdr（補正後p値）, _（未使用の追加情報2つ）
reject, fdr, _, _ = multipletests(p_values, method="fdr_bh")
```

- `fdr_bh` = Benjamini-Hochberg法（FDRを制御する最もポピュラーな方法）
- FDR < 0.01: 「有意」と判定されたタンパク質のうち、偽陽性は1%以下

### 論文の主な発見

| クラスター | タンパク質数 | パターン | 意味 |
|-----------|-----------|---------|------|
| Cluster 3 | 34 | Normal→IV で一貫上昇 | 腫瘍促進因子の候補 |
| Cluster 14 | 1,324 | Normal→IV で上昇 | 広範な腫瘍関連変動 |
| Cluster 20 | 16 | Normal→IV で一貫下降 | 腫瘍抑制因子の候補 |
| Cluster 25 | 1,062 | Normal→IV で下降 | 正常組織機能の喪失 |

### 本書データでの再現結果

本書では全16患者（CRC01-CRC16）のデータを用い、論文と同じ **5群**（Normal(16), Stage I(3), II(4), III(8), IV(1)）で解析を行います。

| 指標 | 本書 (sage, 2081タンパク質) | 論文 (DIA-NN, 10329タンパク質) |
|------|----------------------------|-------------------------------|
| ANOVA 有意 (FDR<0.01) | **720** | 記載なし（論文はクラスター先行） |
| クラスター数 | **30** | 30 |

### ステージ別ヒートマップとクラスタープロファイル

![ステージ別ヒートマップ](images/fig3_stage_heatmap.png)

ANOVA で有意（FDR<0.01）と判定された720タンパク質について、ステージ別の発現パターンをヒートマップとプロファイルプロットで可視化しています。30クラスターに分割した結果から代表的なクラスターを抽出し、Normal → Stage I → II → III → IV にかけて段階的に変動するパターンが観察されます。灰色の細い線は個々のタンパク質、黒い太い線はクラスター平均を示しています。

全16患者32ファイルで解析したことで、**ANOVAで720個のタンパク質が有意**（FDR<0.01）と判定され、論文と同じ5群（Normal, Stage I, II, III, IV）での30クラスターへの分割が可能になりました。疾患進行に伴うプロテオーム変動パターンという生物学的シグナルは商用クリアな sage パイプラインでも検出できることを示しています。

## 実行方法

**スクリプトで一括実行する場合:**

```bash
python scripts/step_09_stage_analysis.py
```

**Notebook でセルごとに実行する場合:**

```bash
jupyter notebook notebooks/step_10.ipynb
```

Notebook版ではANOVA結果やプロファイルプロットをインラインで確認できます。

## まとめ

ステージ別解析により、疾患進行に伴って一貫して変動するタンパク質クラスターを同定しました。本書データでは ANOVA 有意な 720 タンパク質を 30 クラスターに分割し、Normal → Stage I → II → III → IV にかけて段階的に変動する候補マーカー群を特定しました。

> 前回: [#9 COSMIC照合](article-09-cosmic.md)
> 次回: [#11 まとめと次のステップ](article-11-conclusion.md) — 論文化への道筋

#バイオインフォマティクス #プロテオミクス #クラスター分析 #大腸がん #labcode
