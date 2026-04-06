# DIA-MS プロテオミクス × バイオインフォマティクス × 機械学習 バイオマーカー発見論文 10選

> **全てDIA-MS（データ非依存型取得）を起点**とし、「DIA-MS → Spectronaut/DIA-NN → 差分発現 → パスウェイ解析 → 機械学習」のパイプラインでバイオマーカーや薬剤ターゲットを同定した論文をまとめました。X投稿用の要約と使用した解析手法を明記しています。

---

## 技術解説: DIA-MS と 4D-DIA の違い

本まとめでは「DIA-MS」と「4D-DIA」の両方が登場する。以下にその違いを整理する。

### DIA-MS（Data-Independent Acquisition Mass Spectrometry）

質量分析におけるデータ取得戦略の**総称**。従来のDDA（Data-Dependent Acquisition）では強度上位のプリカーサーイオンのみを選択してフラグメンテーションするが、DIA-MSでは設定した質量範囲全体を**網羅的にフラグメンテーション**する。これにより、試料中のペプチドを偏りなく定量できる。

- **分離次元（3次元）:** m/z、保持時間（RT）、フラグメントイオン強度
- **代表的装置:** Thermo Q Exactive シリーズ、Orbitrap Exploris 480、Orbitrap Astral など
- **代表的取得方式:** SWATH-DIA（SCIEX）、通常DIA（Thermo）

### 4D-DIA（イオンモビリティ付きDIA）

DIA-MSに**イオンモビリティ（ion mobility）分離**を第4の次元として追加した手法。Bruker timsTOFシリーズに搭載されたTIMS（Trapped Ion Mobility Spectrometry）技術を利用し、**dia-PASEF**（Parallel Accumulation Serial Fragmentation）方式でデータを取得する。

- **分離次元（4次元）:** m/z、保持時間（RT）、**イオンモビリティ（1/K0）**、フラグメントイオン強度
- **代表的装置:** Bruker **timsTOF Pro / timsTOF HT**
- **代表的取得方式:** **dia-PASEF**

### 比較表

| 特徴 | DIA-MS（従来型） | 4D-DIA |
|------|-----------------|--------|
| 分離次元 | 3次元（m/z, RT, 強度） | **4次元**（m/z, RT, **1/K0**, 強度） |
| 代表的装置 | Q Exactive, Exploris, Astral | Bruker **timsTOF** シリーズ |
| 取得方式 | SWATH-DIA 等 | **dia-PASEF** |
| ピーク分離能 | m/z と RT による分離 | イオンモビリティによる**追加分離**で共溶出ペプチドの識別精度向上 |
| 感度 | 標準的 | TIMSのイオン蓄積効果により**感度向上** |
| 同定数 | 良好 | 一般に**より多くのタンパク質・ペプチドを同定** |
| CCS値 | 取得不可 | **衝突断面積（CCS）**を取得でき、同定の信頼性向上に利用可能 |

### 本まとめにおける使用状況

| 論文 | 装置 | 方式 | タイプ |
|------|------|------|--------|
| 1 (Toyota), 2 (Xing), 8 (Fu), 9 (Yu) | Q Exactive Plus / HF-X | 通常DIA | DIA-MS（3次元） |
| 3 (Gonçalves) | DIA-MS | 通常DIA | DIA-MS（3次元） |
| 4 (Cheng) | DIA-MS（装置詳細なし） | 通常DIA | DIA-MS（3次元） |
| 5 (Martinez-Val) | Orbitrap Exploris 480 + FAIMS | DIA + FAIMS | DIA-MS（3次元 + イオンフィルタリング） |
| 6 (Zhao) | Orbitrap Astral | 通常DIA | DIA-MS（3次元） |
| **7 (Lin)** | **timsTOF（推定）** | **dia-PASEF（推定）** | **4D-DIA（4次元）** |
| **10 (Yin)** | **timsTOF Pro** | **dia-PASEF** | **4D-DIA（4次元）** |

> **ポイント:** 4D-DIAはDIA-MSの発展形であり、イオンモビリティという第4の分離次元を加えることで感度・同定数・定量精度がいずれも向上する。本まとめの10論文中、論文7と10が4D-DIA（timsTOF/dia-PASEF）を使用し、残りの8論文はOrbitrap系の従来型DIA-MSを使用している。なお、論文5のFAIMSもイオンフィルタリングを行うが、TIMSとは原理が異なり、CCS値は取得できないため「4D-DIA」には分類されない。

---

## 1. 【代表論文】大腸がんのがん関連タンパク質同定 (DIA-MS × COSMIC)

**論文:** Toyota N, Konno R, Iwata S, et al. "Identification of Cancer-Associated Proteins in Colorectal Cancer Using Mass Spectrometry." *Proteomes* 2025; 13(3):38. **IF: 3.0**
**DOI:** https://doi.org/10.3390/proteomes13030038

**概要:**
16人のCRC患者から正常/腫瘍組織を採取し、DIA-MSで10,329タンパク質を同定。COSMICデータベースと照合し531のがん関連タンパク質、うち48がCRC特異的。疾患ステージ進行に伴い一貫して増減するタンパク質クラスターを発見。

**パイプライン:** DIA-MS → COSMIC照合 → 差分発現解析 → 階層的クラスタリング → PCA

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| COSMIC | がん関連タンパク質データベース | 無料（学術） |
| R | 統計解析・PCA・クラスタリング | 無料（OSS） |

**X投稿案:**
```
📊 DIA-MSで大腸がんの全プロテオーム解析！

16患者の正常/腫瘍組織から10,329タンパク質を同定
→ COSMICと照合: 531がん関連、48 CRC特異的タンパク質
→ ステージ進行に伴う発現変動クラスターを発見

🔬 解析パイプライン:
DIA-MS(網羅的タンパク質定量)
→ COSMIC(がん体細胞変異DBと照合し、がん関連タンパク質を抽出)
→ 差分発現解析(正常vs腫瘍の発現変動を統計検定)
→ 階層的クラスタリング(発現パターンが類似するタンパク質をグループ化)
→ PCA(高次元データを主成分に圧縮し、サンプル間の全体傾向を可視化)

Toyota et al., Proteomes 2025
https://doi.org/10.3390/proteomes13030038
```

---

## 2. 肝細胞がん早期診断: DIA-MS + Spectronaut → Random Forest → PRM検証

**論文:** Xing X, Cai L, Ouyang J, et al. "Proteomics-driven noninvasive screening of circulating serum protein panels for the early diagnosis of hepatocellular carcinoma." *Nat Commun* 2023; 14:8392. **IF: 14.7**
**DOI:** https://doi.org/10.1038/s41467-023-44255-2

**概要:**
1,002人の血清をDIA-MS(Q Exactive Plus)で探索→Spectronaut(v13.2)でライブラリ構築→Random Forestで4タンパク質パネル(P4)を構築→PRM(Skyline)で検証。HCC vs 肝硬変 AUC 0.979。肝硬変→HCC転換を画像診断の11.4ヶ月前に予測。

**パイプライン:** DIA-MS → Spectronaut → Random Forest → PRM(Skyline)検証
**サンプル:** Discovery 320人 → Retrospective 429人 → Prospective 253人 (計1,002人)

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| Spectronaut v13.2 | DIA-MSデータ解析・ライブラリ構築 | **有料**（Biognosys） |
| Skyline v3.6.0 | PRMデータ解析 | 無料（OSS, MacCoss Lab） |
| R (caret) | Random Forest モデル構築 | 無料（R package） |
| R (pROC) | ROC曲線・AUC算出 | 無料（R package） |

**バイオマーカー:** P4パネル (HABP2, CD163, AFP, PIVKA-II)
- HCC vs LC: AUC 0.979, 感度 0.925
- HCC vs HC: AUC 0.992, 感度 0.975
- LC→HCC転換予測: AUC 0.890（11.4ヶ月前検出）

**X投稿案:**
```
🔬 DIA-MS → Spectronaut → RF で肝がん超早期診断パネル

1,002人の段階的検証:
DIA-MS(探索) → Spectronaut → RF → PRM(Skyline)検証

P4パネル: HABP2, CD163, AFP, PIVKA-II
・HCC vs 肝硬変: AUC 0.979
・HCC vs 健常: AUC 0.992
・肝硬変→HCC: 画像より11.4ヶ月早く予測！

🔬 解析パイプライン:
DIA-MS/Q Exactive Plus(血清タンパク質を網羅的に測定)
→ Spectronaut v13.2(DIAデータからタンパク質を同定・定量しライブラリ構築)
→ Random Forest/caret(決定木アンサンブルで最適マーカーパネルを選定)
→ pROC(ROC曲線・AUCを算出し診断性能を評価)
→ PRM/Skyline(標的タンパク質を高精度に絶対定量し独立コホートで検証)

Xing et al., Nat Commun 2023
https://doi.org/10.1038/s41467-023-44255-2
```

---

## 3. パンキャンサー949細胞株プロテオームマップ: DIA-MS × DIA-NN × 深層学習 → 薬剤脆弱性

**論文:** Gonçalves E, Poulos RC, Cai Z, et al. "Pan-cancer proteomic map of 949 human cell lines." *Cancer Cell* 2022; 40(8):835-849.e8. **IF: 48.8**
**DOI:** https://doi.org/10.1016/j.ccell.2022.06.010
**PMID:** 35839778 | **PMC:** PMC9387775（CC BY オープンアクセス）

**概要:**
949のヒトがん細胞株（28組織、40+がん種）のプロテオームをDIA-MSで6,864ラン測定→DIA-NN(v1.8)でin silicoライブラリを用いて8,498タンパク質を定量。マルチオミクス統合（トランスクリプトーム、薬剤応答625薬剤、CRISPR-Cas9必須遺伝子スクリーン）と深層学習パイプラインにより、転写レベルでは検出できない数千のタンパク質バイオマーカーを同定。タンパク質レベルでのがん脆弱性を体系的に明らかにした大規模リソース（ProCan-DepMapSanger）。

**パイプライン:** DIA-MS → DIA-NN(v1.8, in silicoライブラリ) → マルチオミクス統合(薬剤応答 + CRISPR) → 深層学習 → 薬剤脆弱性予測

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| DIA-NN v1.8 | DIA-MSデータ解析（in silicoライブラリ） | 無料（OSS） |
| Python (深層学習) | 薬剤感受性予測モデル | 無料（OSS） |
| R | 統計解析・マルチオミクス統合 | 無料（OSS） |
| DepMap / GDSC | 薬剤応答・CRISPR必須遺伝子データ | 無料（Web） |
| ProteomeXchange | データ公開 (PXD030304) | 無料（Web） |

**薬剤脆弱性:** 625薬剤 × 947細胞株のスクリーニングと統合し、タンパク質レベルでのみ検出可能ながん脆弱性バイオマーカーを数千個同定

**X投稿案:**
```
🎯 949がん細胞株のパンキャンサーDIA-MSプロテオームマップ

DIA-MS × DIA-NN(v1.8)で6,864ラン→8,498タンパク質を定量
28組織・40+がん種をカバー

625薬剤の応答 + CRISPR必須遺伝子と統合
→ 深層学習で転写レベルでは見えない薬剤脆弱性を発見

全ソフトウェアOSS、データ公開済み(PXD030304)

🔬 解析パイプライン:
DIA-MS(949細胞株を6反復で網羅的タンパク質定量)
→ DIA-NN v1.8(深層学習ベースのOSSでin silicoライブラリからタンパク質同定・定量)
→ マルチオミクス統合(トランスクリプトーム・薬剤応答・CRISPRデータを統合)
→ 深層学習(タンパク質発現パターンからがん薬剤脆弱性を予測)
→ DepMap/GDSC(CRISPR必須遺伝子・薬剤感受性の公開データベースと照合)

Gonçalves et al., Cancer Cell 2022 (CC BY OA)
https://doi.org/10.1016/j.ccell.2022.06.010
```

---

## 4. ベーチェット病: DIA-MS × XGBoost で診断・重症度層別化

**論文:** Cheng L, Li M, Bai Z, et al. "Artificial Intelligence-Driven Proteomics Identifies Plasma Protein Signatures for Diagnosis and Stratification of Behçet's Disease." *Adv Sci* 2025; 12(35):e10061. **IF: 14.3**
**PMID:** 40549882

**概要:**
血漿をDIA-MSで網羅的解析→PLS-DAで可視化→GO/PPIエンリッチメント→XGBoostで診断モデル(AUC 0.984/0.967)→抗体マイクロアレイで検証。重症度層別化(軽度/中等度/重度)にも成功(AUC 0.718-0.960)。補体C4B, C5, C8Aが重症化と関連。初のAIベースBD診断モデル。

**パイプライン:** DIA-MS → PLS-DA → 差分発現 → GO/PPIエンリッチメント → XGBoost → SHAP → 抗体マイクロアレイ検証

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| XGBoost | 診断・層別化モデル | 無料（OSS） |
| SHAP | モデル解釈・特徴量重要度 | 無料（Python package） |
| R / Python | PLS-DA・統計解析 | 無料（OSS） |
| GO / PPI解析ツール | パスウェイ・ネットワーク | 無料（Web） |

**バイオマーカー:** Complement C4B, C5, C8A
- 診断: AUC 0.984 (training), 0.967 (validation)
- 重症度層別化: AUC 0.718-0.960

**X投稿案:**
```
🧫 DIA-MS × XGBoostでベーチェット病の診断・重症度分類に成功

血漿DIA-MSプロテオームにXGBoostを適用:
・診断: AUC 0.984 → 0.967 (検証)
・重症度層別化: AUC 0.718-0.960

補体C4B, C5, C8Aが重症化の鍵
抗体マイクロアレイで独立検証

初のAIベースBD診断・層別化モデル

🔬 解析パイプライン:
DIA-MS(血漿タンパク質を網羅的に定量)
→ PLS-DA(教師あり多変量解析で疾患群vs健常群の分離を可視化)
→ 差分発現解析(BD患者で有意に変動するタンパク質を統計的に同定)
→ GO/PPIエンリッチメント(変動タンパク質の生物学的機能・相互作用を解析)
→ XGBoost(勾配ブースティングで高精度な診断・重症度分類モデルを構築)
→ SHAP(各タンパク質のモデルへの寄与度を定量化し解釈性を付与)
→ 抗体マイクロアレイ(多数の抗体を用いた独立検証で候補マーカーを確認)

Cheng et al., Adv Sci 2025
```

---

## 5. 大腸がん分子サブタイプ: DIA-MS × Spectronaut × limma → CAVIN1再発マーカー

**論文:** Martinez-Val A, et al. "Proteomics of colorectal tumors identifies the role of CAVIN1 in tumor relapse." *Mol Syst Biol* 2025. **IF: 8.3**

**概要:**
371人CRC患者412検体をDIA-MS(Orbitrap Exploris 480 + FAIMS)→Spectronaut v16(directDIA)で解析→ConsensusClusterPlusで4分子サブタイプ同定→limmaで差分発現→ssGSEA/RoKAIでエンリッチメント。EMT様サブタイプでCAVIN1が再発リスクと相関。転写サブタイプ(CMS)より強い予後予測能。siRNAで機能検証。

**パイプライン:** DIA-MS(Exploris 480) → Spectronaut(directDIA) → ComBat → ConsensusClusterPlus → limma → ssGSEA → RoKAI → siRNA検証

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| Spectronaut v16 | DIA-MSデータ解析 (directDIA) | **有料**（Biognosys） |
| Perseus v1.6.15.0 | データフィルタリング | 無料（OSS, MaxQuant系） |
| R v4.1.3 / RStudio | 統計解析全般 | 無料（OSS） |
| ConsensusClusterPlus | 分子サブタイプ分類 | 無料（Bioconductor） |
| limma | 差分発現解析 | 無料（Bioconductor） |
| ComBat (sva) | バッチ補正 | 無料（Bioconductor） |
| GSVA (ssGSEA) | 遺伝子セットエンリッチメント | 無料（Bioconductor） |
| RoKAI | キナーゼ活性推定 | 無料（Web） |

**バイオマーカー:** CAVIN1（再発リスク）, PRELP, ITGA5 / mTOR経路低下が再発と関連

**X投稿案:**
```
🔬 DIA-MS × Spectronaut → limma → ssGSEA で大腸がん再発マーカー

371患者412検体をOrbitrap Exploris 480(DIA)→Spectronaut v16で解析

4分子サブタイプ同定(ConsensusClusterPlus)
→ EMT様サブタイプでCAVIN1が再発リスクと相関
→ CMS転写分類より強い予後予測能
→ siRNAで機能検証済み

🔬 解析パイプライン:
DIA-MS/Exploris 480+FAIMS(イオンフィルタリング付き網羅的タンパク質定量)
→ Spectronaut v16/directDIA(ライブラリ不要でDIAデータからタンパク質を同定・定量)
→ ComBat/sva(バッチ間の技術的偏りを統計補正)
→ ConsensusClusterPlus(複数回のクラスタリングで安定な分子サブタイプを同定)
→ limma(線形モデルで群間の差分発現タンパク質を統計検定)
→ ssGSEA/GSVA(サンプルごとにパスウェイ活性スコアを算出)
→ RoKAI(リン酸化プロテオームからキナーゼ活性を推定)
→ siRNA(遺伝子発現を抑制し候補マーカーの機能を実験的に検証)

Martinez-Val et al., Mol Syst Biol 2025
```

---

## 6. 重症熱性血小板減少症(SFTS): DIA × DIA-NN × XGBoost/RF で予後マーカー

**論文:** Zhao C, Ge Z, Wang R, et al. "Serum proteomics and machine learning identify PSMD11 as a prognostic biomarker in severe fever with thrombocytopenia syndrome." *Front Immunol* 2025; 16:1693946. **IF: 5.7**

**概要:**
SFTS患者87人+健常10人の血清をDIA(Orbitrap Astral)→DIA-NN(v1.8.1, library-free)で解析→642の差分発現タンパク質を同定→GO/KEGG/PPIエンリッチメント→XGBoost+RFで4バイオマーカーを選定。PSMD11がトップ(AUC 0.847)。プロテアソーム経路の異常が致死的転帰と関連。

**パイプライン:** DIA(Orbitrap Astral) → DIA-NN(library-free) → 差分発現(642 DAPs) → GO/KEGG → PPI(Markovクラスタリング) → XGBoost + RF

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| DIA-NN v1.8.1 | DIA-MSデータ解析 (library-free) | 無料（OSS） |
| XGBoost | バイオマーカー選定 | 無料（OSS） |
| Random Forest (R/Python) | バイオマーカー選定 | 無料（OSS） |
| GO / KEGG | パスウェイエンリッチメント | 無料（Web） |
| STRING + Markovクラスタリング | PPIネットワーク解析 | 無料（Web） |

**バイオマーカー:**
- PSMD11: AUC 0.847（LDH, AST, thrombin timeと強相関）
- IL1RL1: AUC 0.847
- PSMC4: AUC 0.843
- IFIH1: AUC 0.791

**X投稿案:**
```
🦠 DIA × DIA-NN × XGBoost/RFでSFTS予後バイオマーカー発見

87患者の血清をOrbitrap Astral(DIA)→DIA-NN(library-free)で解析
→ 642差分発現タンパク質 → GO/KEGG/PPI

XGBoost + RF で4マーカー選定:
PSMD11(AUC 0.847), IL1RL1, PSMC4, IFIH1

プロテアソーム経路の破綻が致死転帰と関連

🔬 解析パイプライン:
DIA/Orbitrap Astral(次世代質量分析計で血清タンパク質を網羅的に定量)
→ DIA-NN v1.8.1(深層学習ベースのOSSでライブラリ不要のタンパク質同定・定量)
→ 差分発現解析(生存群vs死亡群で有意に変動する642タンパク質を同定)
→ GO/KEGG(変動タンパク質の生物学的機能・代謝パスウェイを解析)
→ STRING+Markovクラスタリング(PPIネットワークを構築し機能クラスターに分割)
→ XGBoost(勾配ブースティングで予後予測に重要なマーカーをランク付け)
→ Random Forest(決定木アンサンブルでXGBoost結果をクロスバリデーション)

Zhao et al., Front Immunol 2025
```

---

## 7. ARDS予後予測: 4D-DIA × Spectronaut × LASSO × 5 ML → PRM検証

**論文:** Lin M, Xu F, Sun J, et al. "Integrative multi-omics analysis unravels the host response landscape and reveals a serum protein panel for early prognosis prediction for ARDS." *Crit Care* 2024; 28:223. **IF: 8.8**
**DOI:** https://doi.org/10.1186/s13054-024-05000-3

**概要:**
ARDS患者130人+対照66人の血清を4D-DIA→Spectronaut(v17.5)で2,669タンパク質を定量→KEGG/GSEA→LASSO+Borutaで36候補→5 MLモデル比較→GLMが最良(AUC 0.893)→外部183人でPRM検証し8タンパク質パネルを確立(AUC 0.802)。スフィンゴ脂質シグナリング経路がARDS病態の中心。

**パイプライン:** 4D-DIA → Spectronaut(v17.5) → 差分発現(FC≥2, FDR<0.05) → KEGG/GSEA → PPI(STRING, Cytoscape) → LASSO + Boruta → 5 ML比較(NB, RF, GLM, SVM, GBM) → PRM検証(外部コホート)

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| Spectronaut Pulsar v17.5 | DIA-MSデータ解析 | **有料**（Biognosys） |
| Progenesis QI V2.3 | メタボロミクス解析 | **有料**（Waters/Nonlinear Dynamics） |
| R | LASSO/Boruta/MLモデル全般 | 無料（OSS） |
| glmnet | LASSO回帰 | 無料（R package） |
| Boruta | 特徴量選択 | 無料（R package） |
| STRING | PPIデータベース | 無料（Web） |
| Cytoscape | ネットワーク可視化 | 無料（OSS） |
| KEGG / GO | パスウェイエンリッチメント | 無料（Web） |

**バイオマーカー(8タンパク質パネル):**
- 上昇: VCAM1, LDHB, MSN, FLG2, LMNA, LBP
- 低下: TAGLN2, MBL2
- Discovery AUC: 0.893 / Validation AUC: 0.802

**X投稿案:**
```
🫁 4D-DIA × Spectronaut × 5 ML比較でARDS予後パネル開発

4D-DIA→Spectronaut v17.5で2,669タンパク質定量
→ KEGG/GSEA → LASSO+Boruta → 5 ML比較
→ GLM最良: AUC 0.893 (discovery)

外部183人でPRM検証: AUC 0.802
8タンパク質パネル: VCAM1, LDHB, MSN, FLG2等

スフィンゴ脂質経路がARDS病態の中心

🔬 解析パイプライン:
4D-DIA/timsTOF(イオンモビリティ付きdia-PASEFで血清タンパク質を高感度定量)
→ Spectronaut v17.5(DIAスペクトルからタンパク質を同定・定量)
→ 差分発現解析(FC≥2, FDR<0.05で有意に変動するタンパク質を抽出)
→ KEGG/GSEA(代謝パスウェイ・遺伝子セット単位で機能変動を評価)
→ STRING/Cytoscape(PPIネットワークを構築し中心的ハブタンパク質を特定)
→ LASSO/glmnet(L1正則化で冗長変数を除外しマーカー候補を絞込)
→ Boruta(ランダムフォレストベースの特徴量選択でLASSO結果を補完)
→ 5 ML比較: NB, RF, GLM, SVM, GBM(複数アルゴリズムで最適モデルを選定)
→ PRM/Skyline(標的タンパク質を高精度に定量し外部コホート183人で検証)

Lin et al., Crit Care 2024
```

---

## 8. 胆道閉鎖症早期診断: DIA × Spectronaut × GO/KEGG → Bayesian LR → PRM検証

**論文:** Fu M, Guo Z, Chen Y, et al. "Proteomics Defines Plasma Biomarkers for the Early Diagnosis of Biliary Atresia." *J Proteome Res* 2024; 23:1960-1970. **IF: 3.8**
**PMID:** PMC11077583

**概要:**
新生児(生後60日未満)の血漿をDIA-MS(Q Exactive HF-X)→Spectronaut Pulsar Xで解析→GO/KEGGエンリッチメント→アンサンブル特徴量選択+Bayesianロジスティック回帰で2タンパク質パネル(PIGR, IGL)を構築→PRMで検証。AUC 0.944で胆道閉鎖症を早期鑑別。

**パイプライン:** DIA-MS(Q Exactive HF-X) → Spectronaut Pulsar X → GO/KEGG → PCA → アンサンブル特徴量選択(5-fold CV) → Bayesian LR → PRM検証

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| MaxQuant v1.5.3.17 | DDAデータ解析（補助） | 無料（OSS） |
| Spectronaut Pulsar X | DIA-MSデータ解析 | **有料**（Biognosys） |
| R | Bayesian LR・特徴量選択 | 無料（OSS） |
| GO / KEGG | パスウェイエンリッチメント | 無料（Web） |

**バイオマーカー:** PIGR(上昇) + IGL(低下)
- AUC: 0.944, 感度 0.8, 特異度 1.0

**サンプル:** Discovery 70人(34 BA / 36 IHC) → Validation 84人(49 BA / 35 IHC)

**X投稿案:**
```
👶 DIA × Spectronaut × Bayesian LRで胆道閉鎖症の早期診断

生後60日未満の新生児血漿を
DIA-MS → Spectronaut Pulsar X で解析
→ GO/KEGG → アンサンブル特徴量選択 → Bayesian LR

2タンパク質パネル: PIGR↑ + IGL↓
AUC 0.944、特異度1.0

PRMで独立コホート検証済み
新生児の早期手術介入判断を支援

🔬 解析パイプライン:
DIA-MS/Q Exactive HF-X(新生児血漿タンパク質を網羅的に定量)
→ Spectronaut Pulsar X(DIAデータからタンパク質を同定・定量)
→ MaxQuant(DDAデータを解析しスペクトルライブラリ構築を補助)
→ GO/KEGG(変動タンパク質の生物学的機能・代謝経路を解析)
→ PCA(サンプル間の全体的な発現パターンの違いを可視化)
→ アンサンブル特徴量選択(5-fold CVで複数手法を組合せ安定なマーカーを選定)
→ Bayesian LR(ベイズ推定で事前情報を活用した診断モデルを構築)
→ PRM/Skyline(標的タンパク質を高精度に定量し独立コホートで検証)

Fu et al., J Proteome Res 2024
```

---

## 9. 大腸腺腫の血清診断マーカー: DIA × DIA-NN × LASSO

**論文:** Yu C, Huang X, Cao Y, et al. "Quantitative proteomic analysis reveals potential serum diagnostic markers for colorectal adenoma." *Front Mol Biosci* 2025; 13:1628587. **IF: 3.9**
**DOI:** https://doi.org/10.3389/fmolb.2025.1628587

**概要:**
大腸腺腫患者/炎症性ポリープ患者/健常者の血清をDIA-MS(Q Exactive Plus)→DIA-NN(v1.8.1)で解析→ANOVA(FDR 5%)→sPLS-DA→LASSO回帰でバイオマーカー選定→Metascapeで機能アノテーション。FLNA(フィラミンA)が腺腫で進行性低下し、大腸腺腫の非侵襲的スクリーニングマーカーとして有望(AUC 0.810)。

**パイプライン:** DIA-MS(Q Exactive Plus) → DIA-NN(v1.8.1) → ANOVA(FDR 5%) → sPLS-DA → LASSO → Metascapeエンリッチメント

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| DIA-NN v1.8.1 | DIA-MSデータ解析 | 無料（OSS） |
| R (mixOmics) | sPLS-DA | 無料（R package） |
| R (glmnet) | LASSO回帰 | 無料（R package） |
| Metascape | 機能アノテーション・エンリッチメント | 無料（Web） |

**バイオマーカー:** FLNA (Filamin A)
- 腺腫 vs 健常: AUC 0.810
- 腺腫 vs 炎症性ポリープ: AUC 0.734

**サンプル:** Cohort 1: 138人(プール) → Cohort 2: 42人(個別検証)

**X投稿案:**
```
🔍 DIA × DIA-NN × LASSO で大腸腺腫の血清マーカー発見

血清をDIA-MS → DIA-NN(v1.8.1)で解析
→ ANOVA → sPLS-DA → LASSO

FLNA(フィラミンA)が腺腫で進行性低下
・腺腫 vs 健常: AUC 0.810
・腺腫 vs 炎症性ポリープ: AUC 0.734

大腸内視鏡の補完的スクリーニングに
前がん病変の非侵襲検出を目指す

🔬 解析パイプライン:
DIA-MS/Q Exactive Plus(血清タンパク質を網羅的に定量)
→ DIA-NN v1.8.1(深層学習ベースのOSSでライブラリ不要のタンパク質同定・定量)
→ ANOVA(3群以上の比較で有意に変動するタンパク質を統計検定, FDR 5%)
→ sPLS-DA/mixOmics(スパース教師あり多変量解析で群間差に寄与する変数を抽出)
→ LASSO/glmnet(L1正則化回帰で冗長変数を除外し最小マーカーセットを選定)
→ Metascape(選定マーカーの生物学的機能・パスウェイを統合アノテーション)

Yu et al., Front Mol Biosci 2025
```

---

## 10. 【4D-DIA特集】血清エクソソーム 4D-DIA × Spectronaut × Random Forest で大腸がん診断

**論文:** Yin H, Xie J, Xing S, et al. "Machine learning-based analysis identifies and validates serum exosomal proteomic signatures for the diagnosis of colorectal cancer." *Cell Rep Med* 2024; 5(8):101672. **IF: 11.7**
**PMID:** 39168094

**概要:**
血清エクソソーム(EV)を4D-DIA(timsTOF Pro)で解析→Spectronatuで854タンパク質を同定→OPLS-DA/LASSOnで特徴量選択→5つのMLアルゴリズムを比較しRandom Forestが最適。PF4+AACTの2マーカーでCRC診断AUC 0.960-0.963。912人のELISAで臨床検証。

> **なぜ4D-DIAか？** 本論文はBruker timsTOF Proを用いたdia-PASEF（4D-DIA）でデータを取得している。従来のDIA-MS（3次元: m/z, RT, 強度）に加え、イオンモビリティ（1/K0）を第4の分離次元として利用することで、共溶出するエクソソーム由来ペプチドの分離能が向上し、低存在量タンパク質の検出感度が改善されている。4D-DIAと従来DIA-MSの詳細な違いは冒頭の「技術解説」セクションを参照。

**パイプライン:** 4D-DIA(timsTOF) → Spectronaut → OPLS-DA → LASSO → Random Forest → ELISA検証
**追加解析:** GSEA, STRING PPI, scRNA-seq(Seurat)で細胞起源特定

**使用ソフトウェア（ライセンス）:**
| ソフトウェア | 用途 | ライセンス |
|------------|------|-----------|
| Spectronaut | DIA-MSデータ解析 | **有料**（Biognosys） |
| SIMCA v14.1 | OPLS-DA | **有料**（Sartorius） |
| R v4.2.3 | 統計解析 | 無料（OSS） |
| mlr3 v0.14.1 | ML (RF等) | 無料（R package） |
| glmnet | LASSO回帰 | 無料（R package） |
| IBM SPSS 26.0 | 統計解析 | **有料**（IBM） |
| Cytoscape 3.8.2 | PPIネットワーク可視化 | 無料（OSS） |
| STRING | PPIデータベース | 無料（Web） |
| Seurat v3.1.1 | scRNA-seq解析 | 無料（R package） |
| MSigDB v7.4 | GSEA遺伝子セット | 無料（Broad Institute） |

**バイオマーカー:** PF4 (CXCL4) + AACT → AUC 0.960 (train), 0.963 (test)
**サンプル:** Discovery 37人 → Train 338人 → Test 328人 → External 246人

**X投稿案:**
```
🩸 4D-DIA × Spectronaut × Random Forestで大腸がんリキッドバイオプシー

血清EVをtimsTOF Pro(4D-DIA)で測定
→ Spectronatuで854タンパク質同定
→ OPLS-DA/LASSO → 5 ML比較 → RF最適

PF4 + AACT: AUC 0.960-0.963
CEA/CA19-9を大幅超え、早期CRC(I-II期)でも有効
912人ELISAで検証済み

🔬 解析パイプライン:
4D-DIA/timsTOF Pro(イオンモビリティ付き網羅的タンパク質定量 ※従来DIA-MSとの違いは冒頭の技術解説参照)
→ Spectronaut(DIAスペクトルからタンパク質を同定・定量)
→ OPLS-DA/SIMCA(教師あり多変量解析で群間差を可視化・特徴量抽出)
→ LASSO/glmnet(L1正則化で冗長な変数を除外しマーカー候補を絞込)
→ Random Forest/mlr3(決定木アンサンブルで診断分類モデルを構築)
→ GSEA/MSigDB(遺伝子セット単位で生物学的パスウェイの活性変動を評価)
→ STRING/Cytoscape(タンパク質間相互作用ネットワークを構築・可視化)
→ ELISA(抗体ベースで候補マーカーの血中濃度を大規模臨床検証)

Yin et al., Cell Rep Med 2024
```

---

## 全10論文の共通パイプライン

```
┌─────────────────────────────────────────────────────────────┐
│                    DIA-MS パイプライン                        │
│                                                             │
│  試料  → DIA-MS測定 → Spectronaut/DIA-NN → 差分発現解析     │
│  (血清/    (timsTOF     (ライブラリフリー     (limma/         │
│   血漿/     Orbitrap      or directDIA)       ANOVA/         │
│   組織)     Astral)                           FC+FDR)        │
│                                                             │
│         → パスウェイ解析 → 機械学習 → 検証                    │
│           (GO/KEGG/        (RF/XGBoost/  (PRM/               │
│            GSEA/            LASSO/        ELISA/              │
│            STRING PPI)      Bayesian LR)  マイクロアレイ)      │
└─────────────────────────────────────────────────────────────┘
```

### 装置・ソフトウェア一覧（ライセンス付き）

| カテゴリ | ツール | ライセンス | 使用論文 |
|---------|--------|-----------|---------|
| **MS装置** | timsTOF Pro (4D-DIA) | **有料**（Bruker） | 7, 10 |
| | Orbitrap Exploris 480 | **有料**（Thermo Fisher） | 5 |
| | Orbitrap Astral | **有料**（Thermo Fisher） | 6 |
| | Q Exactive Plus/HF-X | **有料**（Thermo Fisher） | 1, 2, 8, 9 |
| **DIA解析** | Spectronaut | **有料**（Biognosys） | 2, 5, 7, 8, 10 |
| | DIA-NN | 無料（OSS） | 3, 6, 9 |
| **差分発現** | limma | 無料（Bioconductor） | 5 |
| | ANOVA + FDR | 無料（R） | 9 |
| | FC + FDR | 無料（R） | 6, 7 |
| **特徴量選択** | LASSO (glmnet) | 無料（R package） | 7, 9, 10 |
| | OPLS-DA (SIMCA) | **有料**（Sartorius） | 10 |
| | sPLS-DA (mixOmics) | 無料（R package） | 9 |
| | Boruta | 無料（R package） | 7 |
| **ML分類** | Random Forest (caret/mlr3) | 無料（R package） | 2, 6, 10 |
| | XGBoost | 無料（OSS） | 4, 6 |
| | GLM / SVM / GBM | 無料（R package） | 7 |
| | Bayesian LR | 無料（R） | 8 |
| | 深層学習 (Python) | 無料（OSS） | 3 |
| **パスウェイ** | GO / KEGG | 無料（Web） | 4, 6, 7, 8 |
| | GSEA / ssGSEA (GSVA) | 無料（Bioconductor） | 5, 7, 10 |
| | STRING PPI | 無料（Web） | 6, 7, 10 |
| | Metascape | 無料（Web） | 9 |
| **ネットワーク** | Cytoscape | 無料（OSS） | 7, 10 |
| | Perseus | 無料（OSS） | 5 |
| **検証** | PRM + Skyline | 無料（OSS, MacCoss Lab） | 2, 7, 8 |
| | ELISA | - | 10 |
| | 抗体マイクロアレイ | - | 4 |
| | siRNA機能検証 | - | 5 |
| **統計** | R / RStudio | 無料（OSS） | 全論文 |
| | IBM SPSS | **有料**（IBM） | 10 |
| | Progenesis QI | **有料**（Waters） | 7 |
| | MaxQuant | 無料（OSS） | 8 |

### 有料 vs 無料ソフトウェアまとめ

| | 有料 | 無料 |
|---|---|---|
| **DIA解析** | Spectronaut (論文2,5,7,8,10) | DIA-NN (論文3,6,9) |
| **多変量解析** | SIMCA/OPLS-DA (論文10) | mixOmics/sPLS-DA (論文9) |
| **統計** | SPSS (論文10), Progenesis QI (論文7) | R全般 |
| **PRM解析** | - | Skyline (論文2,7,8) |
| **ネットワーク** | - | Cytoscape, STRING, Perseus |
| **ML** | - | XGBoost, RF, LASSO 全て無料 |
| **パスウェイ** | - | GO, KEGG, GSEA, Metascape 全て無料 |

> **ポイント:** DIA解析のSpectronaut以外は、ほぼ全て無料のOSSで構築可能。DIA-NNを使えばDIA解析も無料化でき、パイプライン全体をコストゼロで再現できる（論文6, 9が実例）。

---

## 全OSSパイプライン: 有料ソフトを使わずに再現する方法

本まとめの10論文で使用されている**有料ソフトウェアは全てOSS代替が存在**する。以下に置換マップと、全OSSで構築した場合のパイプラインを示す。

### 有料→OSS置換マップ

| 有料ソフト | 用途 | OSS代替 | 説明 |
|-----------|------|---------|------|
| **Spectronaut** (Biognosys) | DIA-MSデータ解析 | **DIA-NN** (v1.8+) | 深層学習ベースのDIA解析。ライブラリフリーモードで同等以上の性能。論文6,9で使用実績あり。Demichev et al., *Nat Methods* 2020 |
| **SIMCA** (Sartorius) | OPLS-DA多変量解析 | **ropls** (Bioconductor) | OPLS-DAを完全再現するRパッケージ。Thévenot et al., *J Proteome Res* 2015 |
| | | **mixOmics** (R package) | sPLS-DA等の教師あり多変量解析。論文9で使用実績あり |
| **IBM SPSS** (IBM) | 統計解析全般 | **R** (base + tidyverse) | t検定、χ²検定、回帰分析、生存分析等すべてRで再現可能 |
| **Progenesis QI** (Waters) | メタボロミクス解析 | **MS-DIAL** (RIKEN) | LC-MS/MSメタボロミクスのピーク検出・アライメント・同定。Tsugawa et al., *Nat Biotechnol* 2015 |
| | | **XCMS** (Bioconductor) | LC-MSデータの前処理・ピーク検出。Smith et al., *Anal Chem* 2006 |

### 全OSSパイプライン（推奨構成）

```
┌──────────────────────────────────────────────────────────────────┐
│              全OSS DIA-MS バイオマーカー発見パイプライン             │
│                      💰 ライセンス費用: ¥0                        │
│                                                                  │
│  試料 → DIA-MS測定 → DIA-NN → 差分発現解析 → パスウェイ解析       │
│  (血清/    (timsTOF/   (library-  (limma/       (GO/KEGG/         │
│   血漿/    Orbitrap)    free)     DEqMS/        clusterProfiler/  │
│   組織)                           ANOVA)        fgsea)            │
│                                                                  │
│         → 多変量解析 → 特徴量選択 → 機械学習 → 検証               │
│           (ropls/       (glmnet/    (caret/       (Skyline/       │
│            mixOmics)     Boruta)     mlr3/         PRM)           │
│                                      tidymodels)                  │
│                                                                  │
│  可視化: ggplot2 / ComplexHeatmap / Cytoscape                     │
│  ネットワーク: STRING API / igraph / STRINGdb                     │
│  レポート: R Markdown / Quarto                                    │
└──────────────────────────────────────────────────────────────────┘
```

### 各ステップのOSSツール詳細

| ステップ | 推奨OSS | Rパッケージ/ツール | インストール |
|---------|---------|-------------------|-------------|
| **DIA解析** | DIA-NN v1.8+ | CLI / GUI | `github.com/vdemichev/DiaNN` |
| **DDAライブラリ構築** | MaxQuant / FragPipe | CLI | `maxquant.org` / `fragpipe.nesvilab.org` |
| **品質管理** | PTXQC / rawrr | `PTXQC`, `rawrr` | Bioconductor |
| **バッチ補正** | ComBat | `sva::ComBat()` | Bioconductor |
| **差分発現** | limma / DEqMS | `limma`, `DEqMS` | Bioconductor |
| **多変量解析(OPLS-DA)** | ropls | `ropls::opls()` | Bioconductor |
| **多変量解析(sPLS-DA)** | mixOmics | `mixOmics::splsda()` | Bioconductor |
| **特徴量選択(LASSO)** | glmnet | `glmnet::cv.glmnet()` | CRAN |
| **特徴量選択(Boruta)** | Boruta | `Boruta::Boruta()` | CRAN |
| **ML(RF)** | ranger / caret | `ranger`, `caret` | CRAN |
| **ML(XGBoost)** | xgboost | `xgboost::xgb.train()` | CRAN |
| **ML(GLM/SVM/GBM)** | tidymodels | `parsnip`, `workflows` | CRAN |
| **モデル解釈** | SHAP | `SHAPforxgboost`, `DALEX` | CRAN |
| **ROC/AUC** | pROC | `pROC::roc()` | CRAN |
| **パスウェイ(GO/KEGG)** | clusterProfiler | `clusterProfiler::enrichGO()` | Bioconductor |
| **パスウェイ(GSEA)** | fgsea / GSVA | `fgsea::fgsea()`, `GSVA::gsva()` | Bioconductor |
| **PPI解析** | STRINGdb + igraph | `STRINGdb`, `igraph` | Bioconductor / CRAN |
| **ネットワーク可視化** | Cytoscape | `RCy3` (R連携) | OSS |
| **PRM検証** | Skyline | GUI | `skyline.ms` |
| **メタボロミクス** | MS-DIAL / XCMS | `xcms` | RIKEN / Bioconductor |
| **scRNA-seq** | Seurat | `Seurat::FindMarkers()` | CRAN |
| **ヒートマップ** | ComplexHeatmap | `ComplexHeatmap::Heatmap()` | Bioconductor |
| **統計全般** | R + tidyverse | `dplyr`, `ggplot2` 等 | CRAN |

### Rインストールコマンド（一括）

```r
# CRAN パッケージ
install.packages(c(
  "tidyverse", "caret", "ranger", "xgboost", "glmnet", "Boruta",
  "pROC", "DALEX", "SHAPforxgboost", "igraph", "Seurat",
  "survival", "survminer"
))

# Bioconductor パッケージ
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "limma", "DEqMS", "sva",           # 差分発現・バッチ補正
  "ropls", "mixOmics",               # 多変量解析 (OPLS-DA, sPLS-DA)
  "clusterProfiler", "org.Hs.eg.db", # GO/KEGGエンリッチメント
  "fgsea", "GSVA",                    # GSEA
  "ConsensusClusterPlus",             # 分子サブタイプ分類
  "ComplexHeatmap",                   # ヒートマップ
  "STRINGdb",                         # PPI解析
  "RCy3",                             # Cytoscape連携
  "rawrr", "PTXQC",                   # 品質管理
  "xcms"                              # メタボロミクス
))
```

> **まとめ:** 10論文のパイプラインは**Spectronaut → DIA-NN**、**SIMCA → ropls**、**SPSS → R**、**Progenesis QI → MS-DIAL/XCMS** と置換するだけで、全ステップをOSSのみで再現可能。装置（timsTOF, Orbitrap等）はハードウェアのため代替不可だが、ソフトウェア側のライセンス費用は完全にゼロにできる。

---

## 論文の分類（疾患別）

### がん
| 論文 | 疾患 | MS装置 | DIA解析 | ML | AUC |
|------|------|--------|---------|-----|-----|
| 1 (Toyota) | 大腸がん | DIA-MS | - | - | 代表論文 |
| 2 (Xing) | 肝細胞がん | Q Exactive Plus | Spectronaut | RF | 0.979 |
| 3 (Gonçalves) | 40+がん種(949細胞株) | DIA-MS | DIA-NN | 深層学習 | 薬剤脆弱性 |
| 5 (Martinez-Val) | 大腸がん(サブタイプ) | Exploris 480 | Spectronaut | ConsensusClusterPlus | 予後 |
| 9 (Yu) | 大腸腺腫 | Q Exactive Plus | DIA-NN | LASSO | 0.810 |
| **10 (Yin)** | **大腸がん(EV)** | **timsTOF(4D-DIA)** | **Spectronaut** | **RF** | **0.960** |

### 感染症・免疫疾患
| 論文 | 疾患 | MS装置 | DIA解析 | ML | AUC |
|------|------|--------|---------|-----|-----|
| 4 (Cheng) | ベーチェット病 | DIA-MS | - | XGBoost | 0.967 |
| 6 (Zhao) | SFTS | Orbitrap Astral | DIA-NN | XGBoost/RF | 0.847 |

### 重症疾患・小児疾患
| 論文 | 疾患 | MS装置 | DIA解析 | ML | AUC |
|------|------|--------|---------|-----|-----|
| 7 (Lin) | ARDS | 4D-DIA | Spectronaut | 5 ML比較 | 0.893 |
| 8 (Fu) | 胆道閉鎖症 | Q Exactive HF-X | Spectronaut | Bayesian LR | 0.944 |
