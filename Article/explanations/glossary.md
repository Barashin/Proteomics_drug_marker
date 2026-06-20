# DIA-MSプロテオミクス解析 用語集

## あ行

### ANOVA（Analysis of Variance / 分散分析）
**ひとこと**: 3つ以上のグループ間で平均値に差があるかを調べる統計手法

**詳しくは**: t検定が2群比較なのに対し、ANOVAは多群比較（例：Normal vs Stage1 vs Stage2 vs Stage3）を一度に実行できます。プロテオミクスでは疾患ステージ別解析でよく使用されます。

**登場する場所**: step_08_stage_analysis.py, article-10-stage.md

**参考**: Toyota et al. 論文 Methods の "Statistical Analysis"

---

## か行

### COSMIC（Catalogue of Somatic Mutations in Cancer）
**ひとこと**: 世界最大のがん関連遺伝子変異データベース

**詳しくは**: 英国Wellcome Sanger Instituteが運営する、がん組織で発見された遺伝子変異を収集したデータベース。約700個の「がん遺伝子（cancer genes）」を定義しており、発見されたタンパク質ががんとの関連性があるかを調べるために使用されます。

**登場する場所**: step_07_cosmic_analysis.py のセクション1, article-09-cosmic.md

**参考**: https://cancer.sanger.ac.uk/cosmic

---

## さ行

### Sage-proteomics
**ひとこと**: 次世代のオープンソースプロテオミクス検索エンジン

**詳しくは**: Rust言語で開発された高速なDIA解析ツール。MITライセンスのため商用利用に制限がありません。DIA-NNと同等の性能を持ちながら、完全にオープンソースで利用できます。Library-free DIA解析に対応し、理論スペクトルから直接ペプチド同定を行います。

**登場する場所**: step_04_dia_analysis_sage.py 全体

**参考**: 公式リポジトリ https://github.com/lazear/sage

---

## た行

### DIA（Data-Independent Acquisition）
**ひとこと**: 質量分析において、すべての質量範囲を網羅的に測定する手法

**詳しくは**: 従来のDDA（Data-Dependent Acquisition）が強いシグナルから選択的に測定するのに対し、DIAは質量範囲を細かく分割してすべてを測定します。利点は再現性が高く、低発現タンパク質も見逃さないこと。課題はスペクトルが複雑になり、解析により高度な技術が必要になることです。

**登場する場所**: 全スクリプト・全記事で使用される基本概念

**参考**: Toyota et al. 論文 Introduction の "DIA-MS"

---

### t検定（t-test）
**ひとこと**: 2つのグループの平均値に統計的な差があるかを調べる手法

**詳しくは**: プロテオミクスでは「Normal群 vs Tumor群」のようにサンプルを2群に分けて、各タンパク質の発現量に有意差があるかを調べます。Welchのt検定（等分散を仮定しない）を使用するのが一般的。p値 < 0.05で「統計的に有意な差がある」と判定します。

**登場する場所**: step_06_differential_expression.py のセクション2

**参考**: scipyライブラリの stats.ttest_ind 関数

---

## な行

### 欠損値（Missing Value）
**ひとこと**: 測定装置の検出限界以下で「存在するが測定できない」データ

**詳しくは**: プロテオミクスでは低発現タンパク質が「欠損値」として記録されます。これは完全にランダムな欠損ではなく、「発現量が少ない」という生物学的意味を持つMNAR（Missing Not At Random）型の欠損です。単純な削除ではなく、適切な補完手法（downshift法など）が必要です。

**登場する場所**: step_04_preprocess.py のセクション4

**参考**: Perseus software の "Replace missing values from normal distribution"

---

## は行

### PCA（Principal Component Analysis / 主成分分析）
**ひとこと**: 高次元データを低次元に圧縮して、全体のパターンを可視化する手法

**詳しくは**: 数千〜数万のタンパク質データを2次元や3次元に圧縮し、サンプル間の類似性・相違性を散布図で可視化します。第1主成分（PC1）が最も分散の大きい方向、第2主成分（PC2）が2番目に大きい方向を表します。寄与率（%）はその軸がどれだけの情報を保持しているかを示します。

**登場する場所**: step_05_overview_visualization.py のセクション3, article-07-visualization.md

**参考**: scikit-learnライブラリの PCA クラス

---

### PSM（Peptide-Spectrum Match）
**ひとこと**: 測定されたスペクトルとペプチド配列の対応関係

**詳しくは**: 質量分析で得られた個別のMS/MSスペクトルが、どのペプチド（タンパク質断片）由来かを同定した結果です。1つのスペクトルに対して1つのPSMが生成されます。品質はqvalue（False Discovery Rate）で評価され、qvalue < 0.01（FDR < 1%）の同定のみを信頼できるものとして使用します。

**登場する場所**: step_04_dia_analysis_sage.py の出力結果

**参考**: Sage出力の results.sage.tsv ファイル

---

### ProteomeXchange
**ひとこと**: プロテオミクス実験データの国際公開リポジトリ

**詳しくは**: 質量分析による実験データ（RAWファイル）と解析結果を研究者が公開・共有するためのデータベースネットワークです。論文発表時にデータを公開することが多くのジャーナルで義務化されており、再現性研究やメタ解析の基盤となっています。PXDから始まる一意のアクセサリ番号で管理されます。

**登場する場所**: step_02_download_data.py, article-02-data.md

**参考**: Toyota et al. 論文のData Availability, PXD058672

---

## ま行

### mzML
**ひとこと**: 質量分析データの標準オープンフォーマット

**詳しくは**: XMLベースの質量分析データ形式で、装置メーカーに依存しないオープンスタンダードです。ThermoのRAWファイル、Agilentのdファイルなどの装置特有形式をmzMLに変換することで、異なるソフトウェア間でのデータ交換が可能になります。スペクトル情報、保持時間、質量精度などの情報を含みます。

**登場する場所**: step_03_convert.md, step_04_dia_analysis_sage.py の入力

**参考**: Proteomics Standards Initiative (PSI) の仕様

---

## や行

### UniProt
**ひとこと**: 世界最大のタンパク質配列・機能情報データベース

**詳しくは**: Swiss-Prot（手動キュレーション）とTrEMBL（自動アノテーション）で構成される包括的なタンパク質情報リソースです。各生物種の参照プロテオーム（代表的なタンパク質セット）を提供しており、DIA-MS解析では「答え合わせ」用のデータベースとして使用されます。ヒトの場合は約20,000種類のタンパク質配列を収録しています。

**登場する場所**: step_02_download_data.py, step_04_dia_analysis_sage.py

**参考**: https://www.uniprot.org/

---

## ら行

### LFQ（Label-Free Quantification）
**ひとこと**: 化学的標識を使わずにタンパク質の相対定量を行う手法

**詳しくは**: 安定同位体標識（SILAC、iTRAQ等）を使わず、MS1シグナル強度の比較のみでタンパク質の発現量比較を行う方法。コストが安く、サンプル調製が簡単な反面、測定間の変動が大きくなりやすい特徴があります。現在最も広く使用されている定量プロテオミクス手法です。

**登場する場所**: step_04_dia_analysis_sage.py の出力（lfq.tsv）

**参考**: Sage出力ファイルの説明

---

## わ行

### Volcano Plot（火山プロット）
**ひとこと**: 統計的有意性（p値）と生物学的重要性（効果量）を同時に可視化するグラフ

**詳しくは**: 横軸にlog2 fold change（発現変化倍率）、縦軸に-log10(p値)をとった散布図。火山のような形になることからこの名前が付きました。右上・左上の領域にあるタンパク質が「統計的に有意かつ生物学的に重要な変化」を示します。バイオマーカー候補の選定によく使用されます。

**登場する場所**: step_06_differential_expression.py のセクション4, article-08-differential.md

**参考**: matplotlib/seabornでの可視化例

---

## 英数字

### FDR（False Discovery Rate）
**ひとこと**: 間違った発見（偽陽性）の割合

**詳しくは**: 統計的検定において、有意と判定された結果のうち、実際には間違っているものの割合。プロテオミクスでは通常1%未満（FDR < 0.01）に設定され、100個の有意な結果のうち1個以下が偽陽性という意味になります。Benjamini-Hochberg法による多重検定補正でコントロールされます。

**登場する場所**: step_04_dia_analysis_sage.py, step_06_differential_expression.py

**参考**: statsmodelsライブラリの multipletests 関数

---

### Log2 Fold Change
**ひとこと**: 2つの条件間での発現量変化をlog2スケールで表した値

**詳しくは**: Tumor/Normalの比をlog2で表現したもの。log2 FC = 1なら2倍増加、log2 FC = -1なら1/2に減少を意味します。log2変換により、同じ倍率の増加・減少が対称的に表現され（2倍増 = +1, 1/2減 = -1）、統計解析や可視化が容易になります。

**登場する場所**: step_06_differential_expression.py, Volcano plotの横軸

**参考**: numpy.log2 関数による計算

---

### qvalue
**ひとこと**: FDRを考慮した調整済みp値

**詳しくは**: 多重検定補正を行った後のp値。通常のp値は「その検定が間違っている確率」ですが、qvalueは「そのqvalue以下で有意とした全ての検定のうち、間違っているものの割合（FDR）」を表します。プロテオミクスではqvalue < 0.01の基準で有意判定を行います。

**登場する場所**: Sage出力の results.sage.tsv, step_06_differential_expression.py

**参考**: Benjamini-Hochberg法による多重検定補正