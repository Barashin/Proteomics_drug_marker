# step_07_cosmic_analysis.py の解説

## このコードの役割

COSMIC（Catalogue Of Somatic Mutations In Cancer）データベースに登録された既知のがん関連遺伝子と、本研究で同定されたタンパク質を照合し、カバー率を算出するスクリプトです。論文 Section 3.3 に相当する解析を実行します。

## 全体フロー

```
同定タンパク質データ → COSMIC遺伝子リスト照合 → 集合演算 → カバー率計算 → 可視化・CSV保存
```

| | ファイル | 内容 |
|-|---------|------|
| 入力 | `results/preprocessed_data.csv` | 前処理済みタンパク質同定データ |
| 参照 | COSMIC Cancer Gene Census | 既知のがん関連遺伝子データベース |
| 出力 | `results/tables/cosmic_overlap.csv` | 同定されたCOSMIC遺伝子一覧 |
| 出力 | `results/figures/fig_cosmic_coverage.png` | カバー率可視化棒グラフ |

**所要時間の目安**: 約30秒-1分（データベース照合・可視化）

## セクション別解説

### セクション1: COSMIC遺伝子リスト定義

```python
COSMIC_CRC_GENES = [
    # Wntシグナル経路（大腸がんの約80%で異常）
    "APC",       # 腫瘍抑制遺伝子。Wnt経路の負の制御因子
    "TP53",      # 「ゲノムの守護者」。DNA損傷応答の中心
    "KRAS",      # RAS/MAPK経路のがん遺伝子。大腸がんの約40%で変異
    "BRAF",      # MAPK経路のキナーゼ。V600E変異が有名
    "PIK3CA",    # PI3K/AKT経路。約15%で変異
    # ... （全65遺伝子）
]

COSMIC_ALL_GENES = COSMIC_CRC_GENES + [
    # 白血病・リンパ腫・乳がん・肉腫等の追加遺伝子
    "ABL1",      # BCR-ABL融合で有名。イマチニブの標的
    "EGFR",      # 上皮成長因子受容体
    # ... （全389遺伝子）
]
```

**目的**: 科学文献で検証済みのがん関連遺伝子を網羅的にリスト化し、大腸がん特異的遺伝子と全がん種遺伝子の2カテゴリに整理する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `COSMIC_CRC_GENES` | 大腸がんに特異的に関連する65遺伝子のリスト |
| `COSMIC_ALL_GENES` | 全がん種にわたる389遺伝子のリスト |

**引数の意味**:
- CRC特異的: 大腸がんの発症・進展に直接関与するドライバー遺伝子
- 全がん種: CRC遺伝子＋白血病・乳がん・肉腫等の代表的がん遺伝子
- 各遺伝子にコメント付きで機能・変異頻度・関連がん種を記載

**つまずきやすいポイント**:
- COSMICの完全リストは有料なため、論文準拠のサブセットを使用
- 遺伝子名は公式Gene Symbol（HGNC準拠）で統一
- タンパク質名との表記揺れ（例: ERBB2 vs HER2）に注意

### セクション2: 同定タンパク質読み込み

```python
def load_identified_proteins():
    df = pd.read_csv(os.path.join(RESULTS_DIR, "preprocessed_data.csv"), index_col=0)
    return set(df.index.tolist())
```

**目的**: 前処理済みプロテオミクスデータから同定されたタンパク質（遺伝子シンボル）をPythonのset（集合）として読み込む。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `pd.read_csv(index_col=0)` | CSVの1列目（タンパク質名）をインデックスに指定 |
| `set()` | リストを集合オブジェクトに変換（集合演算のため） |
| `df.index.tolist()` | DataFrameのインデックス（行名）をリストとして取得 |

**引数の意味**:
- `index_col=0`: 1列目（Protein列）を行インデックスとして扱う
- `set`: 順序を持たない重複なし要素集合（積集合・差集合演算が高速）

**つまずきやすいポイント**:
- CSVファイルの1行目はヘッダー、1列目はタンパク質名の形式が前提
- 重複遺伝子名があると自動的に除去される（setの特性）

### セクション3: COSMIC照合解析

```python
def cosmic_analysis(identified_proteins):
    cosmic_all = set(COSMIC_ALL_GENES)
    cosmic_crc = set(COSMIC_CRC_GENES)

    # 集合演算による照合
    overlap_all = identified_proteins & cosmic_all    # 積集合
    overlap_crc = identified_proteins & cosmic_crc    # 積集合

    # 差集合（COSMICにあるが同定されなかった遺伝子）
    not_found_all = cosmic_all - identified_proteins
    not_found_crc = cosmic_crc - identified_proteins

    # カバー率の計算
    coverage_all = len(overlap_all) / len(cosmic_all) * 100
    coverage_crc = len(overlap_crc) / len(cosmic_crc) * 100
```

**目的**: Python集合演算を用いてCOSMIC遺伝子リストと同定タンパク質の重複を効率的に計算し、プロテオミクスの網羅性を定量評価する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `&` 演算子 | 積集合（intersection）計算 |
| `-` 演算子 | 差集合（difference）計算 |
| `len()` | 集合の要素数をカウント |

**引数の意味**:
- `overlap_all`: 全がん種COSMIC遺伝子のうち本研究で同定されたもの
- `overlap_crc`: CRC特異的COSMIC遺伝子のうち本研究で同定されたもの
- `not_found_all/crc`: COSMICにあるが本研究で検出されなかった遺伝子
- `coverage`: カバー率（同定数÷COSMIC登録数×100%）

### セクション4: カバー率可視化

```python
def plot_cosmic_summary(result):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 左パネル: 全がん種
    categories = ["COSMIC\n(All Cancer)", "Identified\nin This Study"]
    values = [result["cosmic_all_total"], len(result["overlap_all"])]
    colors = ["#FFB74D", "#4CAF50"]  # オレンジ（COSMIC）、緑（同定済み）

    ax = axes[0]
    bars = ax.bar(categories, values, color=colors, width=0.5, edgecolor="white")
    ax.set_title(f"Cancer-Associated Proteins\n(Coverage: {result['coverage_all']:.1f}%)")

    # 各棒の上部に数値ラベル表示
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                str(val), ha="center", fontweight="bold")
```

**目的**: 2パネル構成の棒グラフで全がん種・CRC特異的のカバー率を直感的に比較可能な形で可視化する。

**主要な関数・ライブラリ**:

| 名前 | 役割 |
|------|------|
| `plt.subplots(1, 2)` | 1行2列のサブプロット作成 |
| `ax.bar()` | 棒グラフ描画 |
| `ax.text()` | 棒の上に数値ラベル配置 |
| `plt.tight_layout()` | サブプロット間の余白自動調整 |

**引数の意味**:
- `figsize=(12, 5)`: 図の横幅12インチ×縦幅5インチ
- `colors=["#FFB74D", "#4CAF50"]`: オレンジ（参照値）・緑（成果値）
- `width=0.5`: 棒の幅（デフォルト0.8より細く）
- `ha="center"`: 水平方向の配置（center=中央揃え）

**つまずきやすいポイント**:
- `bar.get_height() + 5`: 棒の高さにオフセットを追加して数値の位置調整
- カバー率は小数第1位まで表示（`:.1f`）

## 【深掘り】 COSMIC データベースとは？

- **ひとことで**: 世界最大規模の体細胞変異データベースで、がんゲノムの変異を網羅的に収集・整理したリファレンスDB
- **定義**: Catalogue Of Somatic Mutations In Cancer の略。英国 Wellcome Sanger Institute が運営
- **Cancer Gene Census (CGC)**: COSMICの中核で、がんの因果的遺伝子約750個をキュレーション
- **どんなとき使う**: プロテオミクス・ゲノミクス研究で同定した遺伝子の既知がん関連性を検証
- **ライセンス**: 学術利用無料、商用利用は有料契約が必要
- **参考**: https://cancer.sanger.ac.uk/cosmic

## 【深掘り】 集合演算とは？

- **ひとことで**: 数学の集合論に基づく要素間の論理演算で、共通要素・差異要素を効率的に抽出する計算手法
- **積集合 (A & B)**: 両方の集合に含まれる要素（=重複要素）
- **差集合 (A - B)**: AにあるがBにない要素（=Aのみの要素）
- **和集合 (A | B)**: どちらかに含まれる要素（=全要素、重複除去）
- **Python setの利点**: リスト内包表記やループより高速（ハッシュテーブル実装）
- **参考**: Pythonの集合型（set）ドキュメント

## 【深掘り】 カバー率の意味とは？

- **ひとことで**: 既知のがん関連遺伝子のうち、プロテオミクス解析で実際に検出できた割合（%）
- **定義**: カバー率 = (検出されたCOSMIC遺伝子数 / COSMIC総遺伝子数) × 100
- **高いカバー率の意味**: プロテオミクス手法の網羅性が高い、検出感度が優秀
- **論文ベンチマーク**: 全がん種71%、CRC特異的75%が良好な成績の目安
- **技術的制約**: 低発現・膜タンパク質・転写因子は検出困難でカバー率低下の要因
- **参考**: 質量分析プロテオミクスの検出限界と技術的制約

## 【深掘り】 なぜ一部の遺伝子が検出されないのか？

- **ひとことで**: DIA-MSは高感度だが、タンパク質の物理化学的性質・発現量・組織特異性により全遺伝子を検出できない技術的制約がある
- **発現量の問題**: 転写因子（TP53、SOX9等）は極低発現で質量分析の検出限界以下
- **組織特異性**: 造血系遺伝子（ABL1、FLT3）は大腸がん組織では発現しない
- **膜タンパク質**: FGFR、ERBB2等は可溶化困難でペプチド化・検出が困難
- **ペプチドの飛行性**: トリプシン消化で検出に適したペプチドが生成されない場合
- **翻訳後修飾**: リン酸化・糖鎖修飾により参照スペクトルとの照合が失敗

## よくある質問

**Q: 論文値と本実装の値が違うのはなぜ？**
A: COSMICの完全リストは有料のため、本実装では論文記載の代表的遺伝子のサブセットを使用しています。実際の解析では論文補足データ（Table S10-S13）の完全リストを使用します。

**Q: カバー率70%は高いのか低いのか？**
A: プロテオミクス分野では70%以上は非常に優秀な成績です。技術的制約（低発現、膜タンパク質等）を考慮すると、100%達成は現実的ではありません。

**Q: CRC特異的遺伝子のほうがカバー率が高い理由は？**
A: CRC特異的遺伝子は大腸がん組織で実際に発現している可能性が高く、組織特異的でないがん遺伝子（造血系等）が除外されるためです。

**Q: 集合演算を使う利点は？**
A: リストの2重ループ比較（O(n×m)）より、ハッシュテーブルベースの集合演算（O(n+m)）のほうが高速です。数百-数千遺伝子の照合では大幅な性能向上になります。

**Q: 自分のデータで実行したい場合は？**
A: `COSMIC_CRC_GENES`と`COSMIC_ALL_GENES`を自分の対象がん種に合わせて変更し、`preprocessed_data.csv`のパスを変更してください。

## 出力されるファイルの内容

### cosmic_overlap.csv
```
Gene,Type
APC,CRC_Specific
TP53,CRC_Specific
KRAS,CRC_Specific
EGFR,All_Cancer
BCL2,All_Cancer
```

- **Gene**: 本研究で同定されたCOSMIC遺伝子名
- **Type**: CRC_Specific（CRC特異的）またはAll_Cancer（全がん種）
- ソート順: アルファベット順で整列済み

### fig_cosmic_coverage.png
- **左パネル**: 全がん種COSMIC遺伝子のカバー率
- **右パネル**: CRC特異的COSMIC遺伝子のカバー率
- **棒グラフ**: オレンジ（COSMIC登録数）・緑（本研究検出数）
- **数値ラベル**: 各棒の上部に具体的な遺伝子数を表示

## 重要な注意点

- **データベース更新**: COSMICは年1-2回更新されるため、バージョンによって結果が変わる
- **遺伝子名表記**: 公式Gene Symbol（HGNC準拠）で統一されている必要がある
- **ライセンス制約**: 完全なCOSMICリスト利用には学術・商用ライセンス確認が必要
- **解釈の限界**: 検出されない=がんに無関係ではなく、技術的制約の可能性を考慮

## 次のステップ

この解析結果は以下で利用されます：
- Step 8: `step_08_stage_analysis.py`でのステージ別がん関連遺伝子発現解析
- 論文Discussion: プロテオミクス手法の網羅性評価
- バイオマーカー候補タンパク質の既知がん関連性検証