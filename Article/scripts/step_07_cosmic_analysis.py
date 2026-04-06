#!/usr/bin/env python3
# =============================================================================
# Step 5: COSMICデータベースとの照合（Section 3.3 の再現）
# =============================================================================
#
# 【COSMIC（Catalogue Of Somatic Mutations In Cancer）】とは:
#   世界最大規模の体細胞変異データベースであり、英国 Wellcome Sanger Institute
#   が運営している。がんゲノムの体細胞変異（somatic mutation）を網羅的に収集・
#   整理しており、がん研究者にとって不可欠なリファレンスデータベースである。
#   URL: https://cancer.sanger.ac.uk/cosmic
#
# 【Cancer Gene Census（CGC）】とは:
#   COSMICの中核を成すサブデータベースで、体細胞変異によってがんの発症・進展に
#   因果的に関与することが科学的に立証された遺伝子のキュレーション済みリストで
#   ある。各遺伝子には、変異パターン（ドライバー変異 vs パッセンジャー変異）、
#   関連するがん種、がん遺伝子（oncogene）か腫瘍抑制遺伝子（tumor suppressor）
#   かの分類、関連する分子経路などの情報が付与されている。
#   2024年時点で約750遺伝子以上が登録されている。
#
# 【なぜ一部のCOSMIC遺伝子がDIA-MSで検出されないのか】:
#   DIA-MS（Data-Independent Acquisition Mass Spectrometry）は高感度なプロテ
#   オミクス手法であるが、以下の理由により全てのがん関連タンパク質を検出できる
#   わけではない:
#
#   1. 【発現量の問題】: 転写因子（例: TP53, SOX9）やシグナル伝達因子は
#      細胞内での発現量が極めて低く、質量分析の検出限界以下になることがある。
#   2. 【組織特異性】: 一部の遺伝子は特定の組織（例: 造血系のABL1, FLT3）
#      でのみ発現し、大腸がん組織では発現しない場合がある。
#   3. 【タンパク質の物理化学的性質】: 膜貫通タンパク質（例: FGFR, ERBB2）は
#      可溶化が困難で、ペプチド断片化パターンが検出に不利な場合がある。
#   4. 【ペプチドの飛行性（flyability）】: トリプシン消化後に生成されるペプチド
#      がイオン化しにくい、あるいはスペクトルライブラリに登録されていない場合、
#      DIA解析で同定できないことがある。
#   5. 【翻訳後修飾】: 大幅なリン酸化やグリコシル化を受けたタンパク質は、
#      修飾なしの参照スペクトルとの照合が失敗する可能性がある。
#
# 本スクリプトの処理フロー:
#   1. 前処理済みプロテオミクスデータから同定タンパク質リストを読み込む
#   2. COSMIC Cancer Gene Census の遺伝子リストと照合する
#   3. 全がん種およびCRC（大腸がん）特異的な重複遺伝子を抽出する
#   4. カバー率を算出し、棒グラフで可視化する
#   5. 結果をCSVファイルに保存する
#
# 入力: results/preprocessed_data.csv
# 出力: results/tables/cosmic_overlap.csv
#       results/figures/fig_cosmic_coverage.png
#
# 注意: COSMIC APIは学術利用のみ無料。本スクリプトでは論文の補足データ
#       (Table S8-S13) またはCOSMIC公開情報から取得したリストを使用する。
# =============================================================================

# =============================================================================
# ライブラリのインポート
# =============================================================================
import os          # ファイルパス操作・ディレクトリ作成に使用するOS操作モジュール
import numpy as np  # 数値計算ライブラリ（本スクリプトでは直接使用しないが将来の拡張用）
import pandas as pd  # データフレーム操作ライブラリ（CSV読み書き・テーブル操作の中心）
import matplotlib.pyplot as plt  # グラフ描画ライブラリ（棒グラフによるカバー率可視化に使用）

# --- matplotlib_venn の可用性チェック ---
# matplotlib_venn はオプショナル依存ライブラリであり、環境によってはインストール
# されていない場合がある。インストールされていなくてもスクリプトの主要機能（棒グラフ）
# は動作するように、フラグで制御する。
try:
    from matplotlib_venn import venn2  # 【ベン図ライブラリ】の再インポート試行
    HAS_VENN = True  # インポート成功フラグ: ベン図描画機能が利用可能
except ImportError:
    HAS_VENN = False  # インポート失敗フラグ: ベン図描画機能は利用不可

# =============================================================================
# ディレクトリ設定（プロジェクト構成に基づくパス定義）
# =============================================================================
# 【SCRIPT_DIR】: このスクリプト自身が置かれているディレクトリの絶対パス
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# 【PROJECT_DIR】: プロジェクトルート（scriptsの一つ上の階層）
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")

# 【RESULTS_DIR】: 全ての解析結果を格納するディレクトリ
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")

# 【FIG_DIR】: 図（グラフ画像）を格納するサブディレクトリ
FIG_DIR = os.path.join(RESULTS_DIR, "figures")

# 【TABLE_DIR】: 表（CSVファイル）を格納するサブディレクトリ
TABLE_DIR = os.path.join(RESULTS_DIR, "tables")

# 出力先ディレクトリが存在しない場合は自動作成する（exist_ok=Trueで既存時もエラーにしない）
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)


# =============================================================================
# COSMIC Cancer Gene Census の遺伝子リスト定義
# =============================================================================
#
# 【Cancer Gene Census（CGC）】から取得したがん関連遺伝子のリスト。
# 出典: https://cancer.sanger.ac.uk/census
#
# 論文では748タンパク質がCOSMICに登録されており、そのうち531タンパク質
# （約71%）が本研究のDIA-MSプロテオミクスで同定されたと報告されている。
# ここでは論文で報告されたCRC特異的48タンパク質を含む代表的なリストを使用する。
#
# なお、COSMIC CGCの完全なリストはライセンスの関係で再配布できないため、
# 論文の補足データ（Table S10等）に基づく代表的なサブセットを定義している。

# -----------------------------------------------------------------------------
# 【CRC特異的COSMIC遺伝子】
# 大腸がん（CRC: Colorectal Cancer）に特に関連が深いとされる遺伝子群。
# 論文 Table S10 から抽出した65遺伝子のうち代表的なものを収録。
# これらは大腸がんの発生・進展に直接関与するドライバー遺伝子を中心に構成される。
#
# 主なカテゴリ:
#   - Wntシグナル経路: APC, CTNNB1, RNF43, ZNRF3, AXIN2, AMER1
#   - TP53/細胞周期: TP53, RB1, CDK4, CDK8, CCND1, CHEK2
#   - RAS/MAPK経路: KRAS, NRAS, BRAF（大腸がんで最も頻度の高い変異経路）
#   - PI3K/AKT経路: PIK3CA, PTEN
#   - TGF-βシグナル: SMAD4, SMAD2, TGFBR2, ACVR2A
#   - DNA修復/ミスマッチ修復（MMR）: MSH6, MSH2, MLH1, PMS2, MUTYH, POLE, POLD1
#   - 受容体チロシンキナーゼ: ERBB2, ERBB3, MET, IGF1R, FGFR1-3, ALK, RET, ROS1
#   - エピジェネティクス: CREBBP, EP300, KMT2A, KMT2D, ARID1A
#   - その他のがん抑制遺伝子: BRCA1, BRCA2, NF1, VHL, WT1, DCC, PALB2
# -----------------------------------------------------------------------------
COSMIC_CRC_GENES = [
    # --- Wntシグナル経路の遺伝子（大腸がんの約80%で異常が見られる主要経路）---
    "APC",       # 腫瘍抑制遺伝子。Wnt経路の負の制御因子。大腸がんの最初期変異
    "TP53",      # 「ゲノムの守護者」。DNA損傷応答・アポトーシスの中心的制御因子
    "KRAS",      # RAS/MAPK経路のがん遺伝子。大腸がんの約40%で変異
    "BRAF",      # MAPK経路のセリン/スレオニンキナーゼ。V600E変異が有名
    "PIK3CA",    # PI3K/AKT経路の触媒サブユニット。大腸がんの約15%で変異
    "SMAD4",     # TGF-βシグナルの細胞内メディエーター
    "FBXW7",     # E3ユビキチンリガーゼ。MYC, Cyclin E, NOTCHの分解を制御
    # --- RAS/MAPK経路 ---
    "NRAS",      # KRASと同じRASファミリー。大腸がんでは比較的低頻度
    "PTEN",      # PI3K/AKT経路の負の制御因子（ホスファターゼ）
    "CTNNB1",    # β-カテニン。Wnt経路の中心的エフェクター
    "AMER1",     # APC/β-カテニン分解複合体の構成因子（別名: WTX）
    "ARID1A",    # SWI/SNFクロマチンリモデリング複合体のサブユニット
    "SOX9",      # 転写因子。腸管幹細胞の維持に関与
    "TCF7L2",    # Wnt標的遺伝子の転写活性化因子（別名: TCF4）
    # --- TGF-βシグナル経路 ---
    "ACVR2A",    # アクチビン受容体。MSI-H大腸がんで高頻度に変異
    "TGFBR2",    # TGF-β受容体II型。MSI大腸がんで変異頻度が高い
    # --- DNAミスマッチ修復（MMR）遺伝子（リンチ症候群関連）---
    "MSH6",      # ミスマッチ修復タンパク質。MSH2とヘテロダイマーを形成
    "MSH2",      # ミスマッチ修復の中心的タンパク質
    "MLH1",      # ミスマッチ修復タンパク質。プロモーターメチル化で不活性化されることが多い
    "PMS2",      # MLH1とヘテロダイマーを形成するMMRタンパク質
    # --- Wntシグナル関連（追加）---
    "RNF43",     # Wnt受容体Frizzledのユビキチン化を介した負の制御因子
    "ZNRF3",     # RNF43と類似機能を持つE3ユビキチンリガーゼ
    "AXIN2",     # β-カテニン分解複合体の足場タンパク質
    "DCC",       # Deleted in Colorectal Cancer。ネトリン受容体
    "SMAD2",     # TGF-βシグナルのR-SMAD
    # --- DNA損傷応答・修復 ---
    "ATM",       # DNA二本鎖切断の検知・修復シグナルのマスターキナーゼ
    # --- 受容体チロシンキナーゼ（RTK）---
    "ERBB2",     # HER2。分子標的薬トラスツズマブの標的
    "ERBB3",     # HER3。ERBB2とヘテロダイマーを形成しPI3K経路を活性化
    "MET",       # 肝細胞増殖因子（HGF）の受容体
    "IGF1R",     # インスリン様成長因子受容体
    # --- 細胞接着・転写制御 ---
    "CDH1",      # E-カドヘリン。細胞接着分子。EMT（上皮間葉転換）で消失
    "CDK8",      # メディエーター複合体のキナーゼモジュール
    # --- DNA複製忠実性 ---
    "POLE",      # DNAポリメラーゼε。超変異（ultramutation）表現型と関連
    "POLD1",     # DNAポリメラーゼδ。校正（proofreading）活性を持つ
    "MUTYH",     # 塩基除去修復（BER）経路の酸化損傷修復酵素
    # --- その他のシグナル経路 ---
    "BMP4",      # 骨形成タンパク質4。TGF-βスーパーファミリー
    "GNAS",      # Gsα。Gタンパク質共役型受容体シグナルのαサブユニット
    "IDH1",      # イソクエン酸脱水素酵素1。代謝リプログラミングに関与
    "IDH2",      # イソクエン酸脱水素酵素2（ミトコンドリア型）
    "PDGFRA",    # 血小板由来成長因子受容体α
    # --- 受容体チロシンキナーゼ（RTK）追加 ---
    "FGFR1",     # 線維芽細胞増殖因子受容体1
    "FGFR2",     # 線維芽細胞増殖因子受容体2
    "FGFR3",     # 線維芽細胞増殖因子受容体3
    "ALK",       # 未分化リンパ腫キナーゼ。融合遺伝子として重要
    "ROS1",      # 受容体チロシンキナーゼ。ALKと類似した融合パートナー
    # --- JAK/STAT・NOTCH経路 ---
    "RET",       # RET受容体チロシンキナーゼ。甲状腺がん・大腸がんで変異
    "JAK2",      # ヤヌスキナーゼ2。サイトカインシグナル伝達
    "STAT3",     # シグナル伝達兼転写活性化因子3
    "NOTCH1",    # Notchシグナルの受容体1。細胞運命決定に関与
    "NOTCH2",    # Notchシグナルの受容体2
    # --- エピジェネティクス制御因子 ---
    "CREBBP",    # CBP。ヒストンアセチル転移酵素・転写共活性化因子
    "EP300",     # p300。CREBBPと相同なヒストンアセチル転移酵素
    "KMT2A",     # ヒストンH3K4メチル転移酵素（別名: MLL1）
    "KMT2D",     # ヒストンH3K4メチル転移酵素（別名: MLL4）
    "NF1",       # ニューロフィブロミン。RAS-GTPaseの活性化因子（GAP）
    # --- 古典的がん抑制遺伝子 ---
    "VHL",       # フォン・ヒッペル・リンドウ病原因遺伝子。HIF分解を制御
    "WT1",       # ウィルムス腫瘍1。転写因子
    "RB1",       # 網膜芽細胞腫タンパク質。細胞周期G1/Sチェックポイントの制御因子
    "BRCA1",     # 相同組換え修復（HRR）の中心因子。PARP阻害剤の標的
    "BRCA2",     # 相同組換え修復。RAD51のリクルートに必須
    # --- 細胞周期・DNA修復（追加）---
    "CHEK2",     # チェックポイントキナーゼ2。ATMの下流で機能
    "RAD51",     # 相同組換え修復のリコンビナーゼ
    "PALB2",     # BRCA2のパートナー。相同組換え修復に必須
    "CDK4",      # サイクリン依存性キナーゼ4。G1/S期移行を促進
    "CCND1",     # サイクリンD1。CDK4/6の活性化パートナー
]

# -----------------------------------------------------------------------------
# 【全がん種COSMIC遺伝子】
# Cancer Gene Census に登録されている全がん種にわたる遺伝子リスト。
# 実際のCGCには約750遺伝子以上が登録されているが、ライセンスの制約上、
# ここではCRC特異的遺伝子に加え、代表的な汎がん遺伝子を追加したサブセットを定義。
#
# 追加された遺伝子は、白血病・リンパ腫（ABL1, BCL2, BCL6等）、
# 肉腫（EWSR1, FUS等）、乳がん（ESR1, FOXA1等）など
# 多様ながん種に関連するものを含む。
# -----------------------------------------------------------------------------
COSMIC_ALL_GENES = COSMIC_CRC_GENES + [
    # --- チロシンキナーゼ/シグナル伝達（白血病・リンパ腫関連が多い）---
    "ABL1",      # BCR-ABL融合で有名。イマチニブの標的
    "ABL2",      # ABL1ファミリーのキナーゼ
    "ACKR3",     # 非定型ケモカイン受容体3（別名: CXCR7）
    "AFF4",      # 転写伸長因子複合体のサブユニット
    "AKAP9",     # Aキナーゼアンカータンパク質9
    "AKT1",      # PI3K/AKT経路のセリン/スレオニンキナーゼ
    "AKT2",      # AKT1のアイソフォーム。インスリンシグナルにも関与
    # --- 代謝・構造タンパク質 ---
    "ALDH2",     # アルデヒド脱水素酵素2。アセトアルデヒド代謝
    "ANK1",      # アンキリン1。赤血球膜の骨格タンパク質
    "APC2",      # APCファミリーの第2メンバー。Wnt経路の制御
    # --- 核内受容体・キナーゼ ---
    "AR",        # アンドロゲン受容体。前立腺がんの主要ドライバー
    "ARAF",      # RAFキナーゼファミリー。MAPK経路の構成因子
    "ARFRP1",    # ADP-リボシル化因子関連タンパク質1
    # --- RhoGTPase関連 ---
    "ARHGAP26",  # RhoGTPase活性化タンパク質26
    "ARHGEF12",  # Rhoグアニンヌクレオチド交換因子12
    # --- クロマチンリモデリング ---
    "ARID2",     # SWI/SNF複合体のサブユニット（PBAF特異的）
    "ARID5B",    # AT-richインタラクティブドメイン含有タンパク質5B
    "ASXL1",     # ポリコーム関連タンパク質。骨髄系腫瘍で高頻度変異
    "ASXL2",     # ASXL1のパラログ
    # --- DNA損傷応答 ---
    "ATR",       # ATMと相補的なDNA損傷応答キナーゼ（一本鎖切断に応答）
    "ATRX",      # クロマチンリモデラー。ALT（テロメア代替伸長）と関連
    # --- 免疫・抗原提示 ---
    "B2M",       # β2-ミクログロブリン。MHCクラスIの構成因子。免疫逃避と関連
    "BAP1",      # BRCA1関連タンパク質1。ユビキチンカルボキシ末端加水分解酵素
    # --- アポトーシス・転写制御 ---
    "BCL2",      # 抗アポトーシスタンパク質。濾胞性リンパ腫で転座
    "BCL6",      # 転写抑制因子。びまん性大細胞型B細胞リンパ腫のマスター制御因子
    "BCOR",      # BCL6コリプレッサー
    "BCORL1",    # BCORのパラログ
    "BCR",       # BCR-ABL融合遺伝子のパートナー。セリン/スレオニンキナーゼ
    # --- NF-κBシグナル・アポトーシス ---
    "BIRC3",     # cIAP2。アポトーシス阻害タンパク質
    "BLM",       # ブルーム症候群ヘリカーゼ。DNA修復に関与
    "BMPR1A",    # BMP受容体1A。若年性ポリポーシスと関連
    "BTK",       # ブルトンチロシンキナーゼ。B細胞シグナル。イブルチニブの標的
    "BUB1B",     # 紡錘体チェックポイントキナーゼ。染色体安定性に関与
    # --- イオンチャネル・シグナル ---
    "CACNA1D",   # L型カルシウムチャネルαサブユニット
    "CALR",      # カルレティキュリン。骨髄増殖性腫瘍でエクソン9変異
    "CARD11",    # NF-κBシグナルのアダプタータンパク質
    "CARS1",     # システイニルtRNA合成酵素（融合遺伝子パートナー）
    "CASP8",     # カスパーゼ8。外因性アポトーシス経路のイニシエーター
    # --- 転写因子・CBF複合体 ---
    "CBFB",      # コアバインディングファクターβ。急性骨髄性白血病で転座
    "CBL",       # E3ユビキチンリガーゼ。RTKの分解シグナルに関与
    # --- 細胞周期制御 ---
    "CCNB1IP1",  # サイクリンB1相互作用タンパク質1
    "CCND2",     # サイクリンD2。G1/S期移行を促進
    "CCND3",     # サイクリンD3。リンパ球増殖に特に重要
    "CCNE1",     # サイクリンE1。CDK2と結合してG1/S移行を駆動
    # --- 免疫チェックポイント ---
    "CD274",     # PD-L1。免疫チェックポイント分子。免疫療法の標的
    "CD79A",     # BCR（B細胞受容体）シグナルのIgα
    "CD79B",     # BCRシグナルのIgβ
    "CDC73",     # パラフィブロミン。副甲状腺がんで変異
    # --- 細胞接着・キナーゼ ---
    "CDH11",     # カドヘリン11。骨肉腫での融合遺伝子パートナー
    "CDK12",     # RNA Pol II CTDキナーゼ。DNA修復遺伝子の転写に必須
    "CDK6",      # サイクリン依存性キナーゼ6。CDK4と類似機能
    # --- CDKインヒビター（がん抑制因子）---
    "CDKN1A",    # p21。TP53の転写標的。CDK阻害因子
    "CDKN1B",    # p27。細胞周期G1期の停止に関与
    "CDKN2A",    # p16/INK4a。CDK4/6を阻害。多くのがんで欠失
    "CDKN2B",    # p15/INK4b。TGF-β誘導性のCDK阻害因子
    "CDKN2C",    # p18/INK4c。CDK4/6阻害因子
    # --- 転写因子 ---
    "CEBPA",     # CCAAT/エンハンサー結合タンパク質α。骨髄分化に必須
    "CHD4",      # クロマチンリモデラー（NuRD複合体）
    "CHEK1",     # チェックポイントキナーゼ1。ATR下流のDNA損傷応答
    "CIC",       # カピクア転写リプレッサー。RTK/MAPK経路の標的
    "CIITA",     # MHCクラスIIのマスター転写活性化因子
    # --- その他 ---
    "CKS1B",     # CDKサブユニット1B。多発性骨髄腫で増幅
    "CMPK1",     # シチジン一リン酸キナーゼ1
    "COL1A1",    # I型コラーゲンα1鎖。DFSP（皮膚線維肉腫）の融合パートナー
    "CREB1",     # cAMP応答配列結合タンパク質1
    "CREB3L2",   # CREB3-like 2。線維肉腫の融合遺伝子パートナー
    "CRLF2",     # サイトカイン受容体様因子2。ALL（急性リンパ性白血病）で過剰発現
    # --- 受容体チロシンキナーゼ（RTK）---
    "CSF1R",     # コロニー刺激因子1受容体（M-CSF受容体）
    "CSF3R",     # G-CSF受容体。慢性好中球性白血病で変異
    "CUX1",      # Cut-likeホメオボックス1。骨髄系腫瘍の抑制遺伝子
    "CYLD",      # 脱ユビキチン化酵素。NF-κBの負の制御因子
    "DAXX",      # Death domain関連タンパク質。ATRX複合体の構成因子
    "DDR2",      # ジスコイジンドメイン受容体2。コラーゲン受容体チロシンキナーゼ
    "DDX3X",     # DEAD-boxヘリカーゼ3X。RNA代謝・翻訳制御
    # --- RNA関連 ---
    "DICER1",    # RNase III。miRNAの成熟化に必須の酵素
    "DNMT3A",    # DNAメチル転移酵素3A。AML（急性骨髄性白血病）で高頻度変異
    "DROSHA",    # RNase III。miRNA前駆体のプロセシング
    # --- 成長因子受容体 ---
    "EGFR",      # 上皮成長因子受容体。セツキシマブ・パニツムマブの標的
    "EIF4A2",    # 翻訳開始因子4A2
    "ELF4",      # E74様転写因子4
    "ELL",       # RNA Pol II転写伸長因子
    "EP400",     # E1A結合タンパク質p400。クロマチンリモデラー
    # --- Ephrin受容体 ---
    "EPHA3",     # エフリン受容体A3。軸索ガイダンスにも関与
    "EPHA7",     # エフリン受容体A7
    "EPHB1",     # エフリン受容体B1。腸管上皮の区画化に関与
    # --- RTK・転写因子 ---
    "ERBB4",     # HER4。EGFRファミリーの第4メンバー
    "ERG",       # ETSファミリー転写因子。前立腺がんでTMPRSS2-ERG融合
    "ESR1",      # エストロゲン受容体α。乳がんの主要標的
    # --- ETSファミリー ---
    "ETV1",      # ETSファミリー転写因子。GIST（消化管間質腫瘍）で重要
    "ETV4",      # ETSファミリー転写因子。前立腺がん融合遺伝子
    "ETV5",      # ETSファミリー転写因子
    "ETV6",      # ETSファミリー。TEL-AML1融合として白血病で有名
    # --- RNA結合タンパク質 ---
    "EWSR1",     # ユーイング肉腫の融合遺伝子パートナー（EWS-FLI1）
    "EXT1",      # エクソストシン1。ヘパラン硫酸の生合成。骨軟骨腫
    "EXT2",      # エクソストシン2
    # --- エピジェネティクス ---
    "EZH2",      # ポリコーム抑制複合体2（PRC2）の触媒サブユニット。H3K27me3
    "FAM46C",    # 非標準ポリ(A)ポリメラーゼ。多発性骨髄腫で変異
    # --- ファンコニ貧血経路（DNA修復）---
    "FANCA",     # ファンコニ貧血相補群A
    "FANCC",     # ファンコニ貧血相補群C
    "FANCD2",    # ファンコニ貧血相補群D2。ICL（鎖間架橋）修復の中心
    "FANCE",     # ファンコニ貧血相補群E
    "FANCF",     # ファンコニ貧血相補群F
    "FANCG",     # ファンコニ貧血相補群G
    # --- アポトーシス・細胞接着 ---
    "FAS",       # CD95/APO-1。外因性アポトーシスの死受容体
    "FAT1",      # プロトカドヘリン。Wnt/Hippoシグナルの制御
    "FAT4",      # FAT1ファミリー。Hippoシグナル経路
    # --- 成長因子 ---
    "FGF19",     # 線維芽細胞増殖因子19。内分泌型FGF
    "FGF3",      # 線維芽細胞増殖因子3
    "FGF4",      # 線維芽細胞増殖因子4
    # --- がん抑制遺伝子 ---
    "FH",        # フマル酸ヒドラターゼ。腎がん・平滑筋腫と関連
    "FLCN",      # フォリキュリン。Birt-Hogg-Dube症候群の原因遺伝子
    # --- 転写因子・RTK ---
    "FLI1",      # Friendリューケミア組み込み1。ユーイング肉腫の融合パートナー
    "FLT1",      # VEGF受容体1
    "FLT3",      # FMS様チロシンキナーゼ3。AMLで変異・ITD頻出
    "FLT4",      # VEGF受容体3。リンパ管新生に関与
    # --- 転写因子（Forkheadファミリー）---
    "FOXA1",     # Forkhead box A1。前立腺がん・乳がんのパイオニア因子
    "FOXL2",     # Forkhead box L2。顆粒膜細胞腫で体細胞変異
    "FOXO1",     # Forkhead box O1。アポトーシス・代謝制御
    "FOXP1",     # Forkhead box P1。B細胞リンパ腫に関連
    "FUBP1",     # MYCの転写活性化因子
    "FUS",       # FET（FUS/EWSR1/TAF15）ファミリー。多数の肉腫融合遺伝子
    # --- 種々のがん関連遺伝子 ---
    "GAS7",      # 成長停止特異的7。AMLの融合遺伝子パートナー
    "GATA1",     # 赤血球・巨核球分化のマスター転写因子
    "GATA2",     # 造血幹細胞の維持に必須の転写因子
    "GATA3",     # T細胞・乳腺の分化に必須の転写因子
    # --- Gタンパク質シグナル ---
    "GNA11",     # Gタンパク質αサブユニット11。ブドウ膜黒色腫で変異
    "GNA13",     # Gタンパク質αサブユニット13。リンパ腫で変異
    "GNAQ",      # Gタンパク質αサブユニットq。ブドウ膜黒色腫で変異
    # --- その他 ---
    "GPC3",      # グリピカン3。肝細胞がんのバイオマーカー
    "GRM3",      # 代謝型グルタミン酸受容体3。メラノーマで変異
    # --- ヒストンバリアント ---
    "H3-3A",     # ヒストンH3.3バリアント。小児びまん性内在性橋膠腫で変異
    "H3-3B",     # ヒストンH3.3バリアント（第2の遺伝子座）
    "H3C2",      # 標準型ヒストンH3
    # --- 成長因子 ---
    "HGF",       # 肝細胞増殖因子。METのリガンド
    "HIST1H3B",  # ヒストンH3.1。小児高悪性度神経膠腫で変異
]


# =============================================================================
# 関数定義
# =============================================================================

def load_identified_proteins():
    """
    前処理済みプロテオミクスデータから同定されたタンパク質リストを読み込む。

    【処理内容】:
        results/preprocessed_data.csv を読み込み、インデックス列（遺伝子名）を
        Python の set（集合）オブジェクトとして返す。
        set を使用するのは、後続のCOSMIC照合で集合演算（積集合・差集合）を
        効率的に行うためである。

    【入力ファイル】:
        results/preprocessed_data.csv
        - インデックス列: 遺伝子シンボル（例: "APC", "TP53"）
        - データ列: 各サンプルにおけるタンパク質発現量

    Returns:
        set: 同定されたタンパク質（遺伝子シンボル）の集合
    """
    # CSVファイルを読み込み、最初の列をインデックス（遺伝子名）として使用
    df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)
    # インデックスをリスト化し、setに変換して返す（集合演算のため）
    return set(df.index.tolist())


def cosmic_analysis(identified_proteins):
    """
    COSMIC Cancer Gene Census データベースとプロテオミクス同定タンパク質を照合する。

    【処理内容】:
        1. COSMIC遺伝子リスト（全がん種・CRC特異的）をset化する
        2. 同定タンパク質との【積集合（intersection）】を計算し、重複遺伝子を抽出
        3. COSMIC登録だが未同定の遺伝子（【差集合】）も記録する
        4. 【カバー率】（= 同定された割合）を百分率で算出する

    【背景】:
        このカバー率は、DIA-MSプロテオミクスの網羅性を評価する重要な指標である。
        論文では、全がん種で748中531（71%）、CRC特異的で64中48（75%）の
        カバー率が報告されており、DIA-MSの高い検出能力を示している。

    Args:
        identified_proteins (set): 同定されたタンパク質の遺伝子シンボル集合

    Returns:
        dict: 以下のキーを含む辞書:
            - cosmic_all_total (int): 全がん種COSMIC遺伝子の総数
            - cosmic_crc_total (int): CRC特異的COSMIC遺伝子の総数
            - overlap_all (set): 全がん種で同定された遺伝子の集合
            - overlap_crc (set): CRC特異的で同定された遺伝子の集合
            - not_found_all (set): 全がん種で未同定の遺伝子の集合
            - not_found_crc (set): CRC特異的で未同定の遺伝子の集合
            - coverage_all (float): 全がん種のカバー率（%）
            - coverage_crc (float): CRC特異的のカバー率（%）
    """
    # リストをsetに変換（集合演算を可能にするため）
    cosmic_all = set(COSMIC_ALL_GENES)  # 全がん種のCOSMIC遺伝子セット
    cosmic_crc = set(COSMIC_CRC_GENES)  # CRC特異的のCOSMIC遺伝子セット

    # 【集合演算による照合】
    # & 演算子: 積集合（intersection）= 両方に含まれる遺伝子
    overlap_all = identified_proteins & cosmic_all    # 全がん種で同定された遺伝子
    overlap_crc = identified_proteins & cosmic_crc    # CRC特異的で同定された遺伝子

    # - 演算子: 差集合（difference）= COSMICにあるが同定されなかった遺伝子
    not_found_all = cosmic_all - identified_proteins  # 全がん種で未検出の遺伝子
    not_found_crc = cosmic_crc - identified_proteins  # CRC特異的で未検出の遺伝子

    # 結果を辞書として返す
    return {
        "cosmic_all_total": len(cosmic_all),        # COSMIC全がん種の登録遺伝子数
        "cosmic_crc_total": len(cosmic_crc),        # COSMIC CRC特異的遺伝子数
        "overlap_all": overlap_all,                  # 全がん種で同定された遺伝子セット
        "overlap_crc": overlap_crc,                  # CRC特異的で同定された遺伝子セット
        "not_found_all": not_found_all,              # 全がん種で未検出の遺伝子セット
        "not_found_crc": not_found_crc,              # CRC特異的で未検出の遺伝子セット
        # 【カバー率の計算】: 同定数 / COSMIC登録数 * 100（%）
        # ゼロ除算防止のため、cosmic_allが空の場合は0を返す
        "coverage_all": len(overlap_all) / len(cosmic_all) * 100 if cosmic_all else 0,
        "coverage_crc": len(overlap_crc) / len(cosmic_crc) * 100 if cosmic_crc else 0,
    }


def plot_cosmic_summary(result):
    """
    COSMICデータベースとの照合結果を棒グラフで可視化する。

    【処理内容】:
        2パネル構成の棒グラフを作成する:
        (a) 左パネル: 全がん種のCOSMIC登録数 vs 同定数（カバー率表示付き）
        (b) 右パネル: CRC特異的のCOSMIC登録数 vs 同定数（カバー率表示付き）

    【可視化のポイント】:
        - オレンジ色: COSMICデータベースの登録数（基準値）
        - 緑色: 本研究で同定された数（成果指標）
        - 各棒の上部に具体的な数値を表示
        - タイトルにカバー率を併記

    Args:
        result (dict): cosmic_analysis() の戻り値辞書

    Returns:
        None（図をファイルに保存する副作用のみ）
    """
    # 1行2列のサブプロット構成で図を作成（横12インチ x 縦5インチ）
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # =====================================================================
    # (a) 左パネル: 全がん種のCOSMICカバー率
    # =====================================================================
    categories = ["COSMIC\n(All Cancer)", "Identified\nin This Study"]  # X軸ラベル
    values = [result["cosmic_all_total"], len(result["overlap_all"])]   # 棒グラフの値
    colors = ["#FFB74D", "#4CAF50"]  # オレンジ（COSMIC）、緑（同定済み）

    ax = axes[0]  # 左パネルのAxesオブジェクトを取得
    # 棒グラフの描画（width=0.5で棒幅を調整、edgecolor="white"で視認性向上）
    bars = ax.bar(categories, values, color=colors, width=0.5, edgecolor="white")
    # タイトルにカバー率を含めて設定
    ax.set_title(f"Cancer-Associated Proteins\n(Coverage: {result['coverage_all']:.1f}%)")
    ax.set_ylabel("Number of Proteins")  # Y軸ラベル
    # 各棒の上部に数値ラベルを追加（棒の高さ + 5の位置に太字で表示）
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                str(val), ha="center", fontweight="bold")

    # =====================================================================
    # (b) 右パネル: CRC特異的のCOSMICカバー率
    # =====================================================================
    categories = ["COSMIC\n(CRC)", "Identified\nin This Study"]         # X軸ラベル
    values = [result["cosmic_crc_total"], len(result["overlap_crc"])]   # 棒グラフの値

    ax = axes[1]  # 右パネルのAxesオブジェクトを取得
    # 棒グラフの描画（全がん種と同じ色スキームを使用）
    bars = ax.bar(categories, values, color=colors, width=0.5, edgecolor="white")
    # タイトルにカバー率を含めて設定
    ax.set_title(f"CRC-Associated Proteins\n(Coverage: {result['coverage_crc']:.1f}%)")
    ax.set_ylabel("Number of Proteins")  # Y軸ラベル
    # 各棒の上部に数値ラベルを追加（CRCは値が小さいためオフセットは+2）
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                str(val), ha="center", fontweight="bold")

    # =====================================================================
    # 共通のスタイル設定（両パネルに適用）
    # =====================================================================
    for ax in axes:
        # 上辺と右辺の枠線を非表示にする（学術論文スタイル）
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # サブプロット間の余白を自動調整
    plt.tight_layout()

    # 図をPNGファイルとして保存（150dpi、余白を最小化）
    filepath = os.path.join(FIG_DIR, "fig_cosmic_coverage.png")
    fig.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close()  # メモリ解放のために図を閉じる
    print(f"    保存: {filepath}")


# =============================================================================
# メイン処理
# =============================================================================

def main():
    """
    COSMIC照合解析のメインエントリーポイント。

    【処理フロー】:
        1. 前処理済みプロテオミクスデータから同定タンパク質リストを読み込む
        2. COSMIC Cancer Gene Census との照合を実行し、カバー率を算出する
        3. 照合結果をCSVファイルに保存する
        4. 棒グラフによる可視化を行い、PNG画像として保存する

    【出力ファイル】:
        - results/tables/cosmic_overlap.csv : 同定されたCOSMIC遺伝子の一覧
        - results/figures/fig_cosmic_coverage.png : カバー率の棒グラフ

    【論文との対応】:
        本関数の出力は、論文 Section 3.3（COSMIC Analysis）の Fig.5 および
        Table S10-S11 に対応する。
    """
    # --- ヘッダー表示 ---
    print("=" * 60)
    print("Step 5: COSMICデータベースとの照合")
    print("=" * 60)

    # =====================================================================
    # 1. 同定タンパク質の読み込み
    # =====================================================================
    print("\n[1/3] 同定タンパク質の読み込み")
    identified = load_identified_proteins()  # setとして読み込み
    print(f"  同定タンパク質数: {len(identified)}")

    # =====================================================================
    # 2. COSMIC照合の実行
    # =====================================================================
    print("\n[2/3] COSMIC照合")
    result = cosmic_analysis(identified)  # 集合演算による照合

    # --- 全がん関連タンパク質の結果表示 ---
    print(f"\n  全がん関連タンパク質:")
    print(f"    COSMIC登録: {result['cosmic_all_total']}")          # CGC登録総数
    print(f"    本研究で同定: {len(result['overlap_all'])}")         # 同定された数
    print(f"    カバー率: {result['coverage_all']:.1f}%")            # カバー率（%）
    print(f"    （論文値: 748中531 = 71%）")                         # 論文との比較参考値

    # --- CRC関連タンパク質の結果表示 ---
    print(f"\n  CRC関連タンパク質:")
    print(f"    COSMIC登録: {result['cosmic_crc_total']}")          # CRC CGC登録総数
    print(f"    本研究で同定: {len(result['overlap_crc'])}")         # 同定された数
    print(f"    カバー率: {result['coverage_crc']:.1f}%")            # カバー率（%）
    print(f"    （論文値: 64中48 = 75%）")                           # 論文との比較参考値

    # =====================================================================
    # 結果をCSVファイルに保存
    # =====================================================================
    # 全がん種で重複した遺伝子のDataFrameを作成（アルファベット順にソート）
    overlap_df = pd.DataFrame({
        "Gene": sorted(result["overlap_all"]),                       # 遺伝子名（ソート済み）
        "Type": ["All_Cancer"] * len(result["overlap_all"]),         # 分類ラベル
    })
    # CRC特異的で重複した遺伝子のDataFrameを作成
    crc_df = pd.DataFrame({
        "Gene": sorted(result["overlap_crc"]),                       # 遺伝子名（ソート済み）
        "Type": ["CRC_Specific"] * len(result["overlap_crc"]),       # 分類ラベル
    })
    # 2つのDataFrameを縦方向に結合してCSVに保存
    pd.concat([overlap_df, crc_df]).to_csv(
        os.path.join(TABLE_DIR, "cosmic_overlap.csv"), index=False   # インデックス列は不要
    )

    # =====================================================================
    # 3. 可視化（棒グラフの作成・保存）
    # =====================================================================
    print("\n[3/3] 可視化")
    plot_cosmic_summary(result)  # 棒グラフを作成してPNG保存

    # --- フッター表示 ---
    print("\n" + "=" * 60)
    print("COSMIC解析完了！")
    print("=" * 60)


# =============================================================================
# スクリプト実行のエントリーポイント
# =============================================================================
# このファイルが直接実行された場合（python step_05_cosmic_analysis.py）のみ
# main() を呼び出す。他のスクリプトからimportされた場合は実行しない。
if __name__ == "__main__":
    main()
