---
title: "📰 同様の流れで作成されている論文たち"
free: false
---

以下にターゲットはそれぞれユニークですが、同様の手順で、publishされている論文を紹介します。ここ5年の論文が多いので、まだまだこのやり方で論文は書けそうな感じはしています。一部有料ツールを使っています。
状況に応じて、各論文のMaterial and Methodsを参考にすると良いです。

## 1. **Structural Dynamics of P-Rex1 Complexed with Natural Leads Establishes the Protein as an Attractive Target for Therapeutics to Suppress Cancer Metastasis**

https://onlinelibrary.wiley.com/doi/10.1155/2023/3882081

**出版年**:2023年

**標的**：P-Rex1 (PDB: 5D27)
**標的選定理由**：ガンでの様々な機能
**ライブラリ**：CMNPD, MPD3, selleckchem bioactive library I, selleckchem bioactive library II, Asinex target oncology library
**スクリーニング**：FAF-Drugs4
**分子ドッキング**：PyRx, GOLD
**MD simulation**: AMBER20
**評価**：RMSD, RMSF, RoG, βFactor, MM/GBSA, MM/PBSA, RDF plots,
**ADMET analysis**:Swiss ADME, pkCSM, PreADMET


## 2. **Determination of Novel Anti-Cancer Agents by Targeting OGG1 Enzyme Using Integrated Bioinformatics Methods**

https://www.mdpi.com/1660-4601/18/24/13290

**標的**：OGG1 enzyme (PDB: 6RLW)
**標的選定理由**：OGG1の機能低下が様々な病態を阻害
**ライブラリ**：CMNPD, MDP3, selleckchem bioactive library I, II, Asinex target oncology library
**スクリーニング1**：Discovery Studio 3.5 toxicity prediction module
**スクリーニング2**：PyRx, GOLD
**MD**: AMBER20
**評価**：RMSD, RMSF, RoG, βFactor, MM/GBSA, MM/PBSA, Hydrogen bonds analysis, WaterSwap, Alanine scan
**ADMET analysis**: Swiss ADME, pkCSM


## 3. **In silico screening of potential β-secretase (BACE1) inhibitors from VIETHERB database**

https://doi.org/10.1007/s00894-022-05051-9

**出版年**: 2022年
**IF**: 2.2
**標的**: β-secretase (BACE1) (PDB: 2P4J)
**標的選定理由**: β-amyloid形成に関係、阻害剤はアルツハイマーに効く
**スクリーニング**: AutoDock4
**MD**: GROMACS 5.1.3
**評価法**: Binding free energy, Steered Molecular Dynamics (SMD) simulations, RMSD,
**ADMET**: PreADMET

## 4.  **In Silico Identification of Angiotensin-Converting Enzyme Inhibitory Peptides from MRJP1**

https://doi.org/10.1371/journal.pone.0228265

**出版年**: 2020年
**IF**: 3.24
**標的**: Angiotensin-Converting Enzyme (ACE) (PDB: 1O8A)
**標的選定理由**: ACEは血圧の調節に重要な役割を果たし、その阻害は高血圧を予防するための主要なターゲットとされています。
**ライブラリ**: AHTpin serverを使ったMRJP1由来のペプチドライブラリ
**スクリーニング**: PatchDock
**分子ドッキング**: FireDock
**MD**: AMBER v14
**評価法**: RMSD, RMSF, B-factor, QSARモデル, SVMスコア

## 5. **Molecular Insights into Binding Mode and Interactions of Structure-Based Virtually Screened Inhibitors for *Pseudomonas aeruginosa* Multiple Virulence Factor Regulator (MvfR)**

https://doi.org/10.3390/molecules26226811

**Publication Year**: 2021
**IF**: 4.2
**標的**: Multiple Virulence Factor Regulator (MvfR) (PDB: 6B8A)
**標的選定理由**: MvfR を阻害することで、緑膿菌を抗生剤に効きやすくする。
**ライブラリ:** Asinex Antibacterial Library , CMNPD
**スクリーニング**: AutoDock Vina, GOLD 5.2
**MD**: AMBER20
**評価法**: RMSD, RMSF, RoG, Hydrogen bond analysis, RDF, MM/PBSA and MM/GBSA binding free energies
**ADMET Tools**: SwissADME, pkCSM

## 6. **(単著！！！）Identification of promising methionine aminopeptidase enzyme inhibitors: A combined study of comprehensive virtual screening and dynamics simulation study**

https://doi.org/10.1016/j.jsps.2023.101745

**出版年**: 2023年
**IF**: 3.0
**標的**: Methionine Aminopeptidase (MetAP) enzyme from *Rickettsia prowazekii* (PDB: 3MR1)
**標的選定理由**: 細菌の生存に不可欠、これが抗菌標的として非常に魅力的です。
**ライブラリ**: Asinex Antibacterial Library, CMNPD
**スクリーニング**: PyRx
**MD**: AMBER 22
**評価法**: RMSD, RMSF, Radius of Gyration (Rg), Hydrogen Bond Analysis, MM/GBSA, MM/PBSA, WaterSwap Method
**ADMET**: ADMETlab 2.0

## 7. ***In silico* high-throughput virtual screening and molecular dynamics simulation study to identify inhibitor for AdeABC efflux pump of *Acinetobacter baumannii***

https://www.tandfonline.com/doi/full/10.1080/07391102.2017.1317025

**出版年**: 2017年
**標的**: AdeB protein of AdeABC efflux pump (PDB: Homology model, no PDB ID provided)
**標的選定理由**: AdeABC effluxポンプは、薬剤排出機構によりAcinetobacter baumanniiの薬剤耐性を引き起こすため、薬剤のターゲットとして有効である。
**タンパク質の構造評価 (全無料）**：PROCHECK([https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/](https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/))、VERIFY3D([https://www.doe-mbi.ucla.edu/verify3d/](https://www.doe-mbi.ucla.edu/verify3d/))、Ramachandran  plots、ProSA-Web（[https://prosa.services.came.sbg.ac.at/prosa.php](https://prosa.services.came.sbg.ac.at/prosa.php)）、PDBsum([https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/](https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/)）
**活性ポケット発見**：SiteMap
**ライブラリ(無料)**: ZINC database (BioBlocks library, 159,868 compounds)
**スクリーニング**: Maestro v10.3, QikProp (ADMET analysis), Glide (HTVS, SP, XP docking)
**分子ドッキング**: Glide (Schrödinger suite)
**MD**: GROMACS v5.1.4
**ADMET**: QikProp

## 8. ***In silico* screening of *Andrographis paniculata* secondary metabolites as anti-diabetes mellitus through PDE9 inhibition**

https://journals.lww.com/rips/fulltext/2023/18010/in_silico_screening_of_andrographis_paniculata.10.aspx

**出版年:** 2023年
**IF**: 2.6
**標的**: PDE9 (PDB: 4Y87)
**標的選定理由**: PDE9はcGMPの分解に関与し、インスリン抵抗性の改善に寄与する可能性があるため。
**ライブラリ**: Maebashi Institute of TechnologyおよびNara Institute of Science and Technologyによる*Andrographis paniculata*二次代謝産物ライブラリ　（[http://www.knapsackfamily.com/KNApSAcK_Family/](http://www.knapsackfamily.com/KNApSAcK_Family/)）
**スクリーニング**: AutoDockTools, AutoDock
**MD**: Gromacs 2021.3, AMBER99SB-ILDN force field
**評価法**: RMSD, MM/GBSA, MM/PBSA, Hydrogen bonds analysis

## 9. **MOLECULAR DOCKING AND MOLECULAR DYNAMIC  STUDIES: SCREENING PHYTOCHEMICALS OF Acalypha  indica AGAINST BRAF KINASE RECEPTORS FOR POTENTIAL USE IN MELANOCYTIC TUMOURS**

https://rasayanjournal.co.in/admin/php/upload/3611_pdf.pdf

**IF**: 0.81
**標的**: BRAF kinase receptor (PDB: 6XFP)
**標的選定理由**: BRAFキナーゼは細胞アポトーシスに関与しており、悪性黒色腫における主要な標的となるため。
**ライブラリ**: *Acalypha indica*の植物化学成分（the NCBI PubChem Compound databaseから抽出）
**スクリーニング**: AutoDock 4.2
**MD**: GROMACS 2019.6
**評価法**: RMSD, MM/PBSA法による結合エネルギー、残基の結合エネルギーへの寄与評価

## 10. **Bulbine frutescens phytochemicals as novel ABC-transporter inhibitor: a molecular docking and molecular dynamics simulation study**

http://dx.doi.org/10.20517/2394-4722.2020.92

**出版年**: 2021年
**IF**: 1.9
**標的**: P-glycoprotein (PDB: 3G60)
**標的選定理由**: ABCトランスポーターは、薬剤排出を介して薬剤耐性に関与し、癌治療の効果を低減することが知られているため、ABCトランスポーターの阻害が多剤耐性克服の戦略となり得る。
**ライブラリ**: Bulbine frutescensの既知のフィトケミカルライブラリ 25 個
**スクリーニング**: AutoDock Tools 1.5.6
**MD**: GROMACS (version 2018)
**評価法**: RMSD, RMSF, Radius of Gyration (RoG), 水素結合解析、MM-GBSAによるタンパク質・リガンド結合自由エネルギー計算

## 11. **Identification of lead inhibitors for 3CLpro of SARS-CoV-2 target using machine learning based virtual screening, ADMET analysis, molecular docking and molecular dynamics simulations**

https://pubs.rsc.org/en/content/articlelanding/2024/ra/d4ra04502e

**出版年**: 2024年
**IF**: 4.03
**標的**: 3CLpro (PDB: 5R82)
**標的選定理由**: 3CLproはSARS-CoV-2ウイルスの複製に関与し、ウイルス感染を抑制するための主要なターゲット酵素です。
**ライブラリ**: eMolecules databases https://search.emolecules.com/
**スクリーニング**: CHEMBLからさせたモデルによるバーチャルスクリーニング
**分子ドッキング**: AutoDock 4
**MDシミュレーション**: GROMACS (CHARMM36力場)
**評価法**: RMSD, RMSF, Protein Solvent Accessible Surface Area (SASA), RoG, MM-PBSA (結合自由エネルギー)など
**ADMET**: SwissADME, ProTox-II, pkCSM

## 12. **Discovery of human coronaviruses pan-papain-like protease inhibitors using computational approaches**

https://www.sciencedirect.com/science/article/pii/S2095177920310558?via%3Dihub

**標的**: Papain-like protease (PLpro) (PDB: 6W9C, 3MJ5, 4R3D)
**標的選定理由**: PLproはコロナウイルスの複製と免疫回避に重要な役割を果たすため、抗ウイルス薬のターゲットとして有望
**ライブラリ**: Asinex protease inhibitor library (~7,000化合物)
**スクリーニング**: PyRx 0.8（6W9Cを使う）
**分子ドッキング**:Autodock-Vina 1.1.2
**MD**: GROMACS 2018
**パラメーター**: SwissParamを使って、トポロジーファイル生成、OPLS-AA/L force-field、TIP3P
**評価法**: RMSD, RMSF, ポテンシャルエネルギー, RoG, 水素結合解析, SASA, PCA, MM-PBSA, MM-GBSA
**ADMET**: SwissADME, PreADMET, pkCSM

## 13. **Virtual screening of knottin and defensin peptides perceives hits against the SARS CoV-2 RBD domain and hACE2 interaction**

https://chemrxiv.org/engage/chemrxiv/article-details/63c7a4da23c13b2e4813aa5b

**出版年**: 2023年
**IF**: なし　ChemRxiv
**標的**: SARS-CoV-2 RBD (PDB: 記載なし)
**標的選定理由**: SARS-CoV-2 RBDドメインとヒトACE2の相互作用を阻害し、ウイルスの宿主細胞への侵入を防ぐことを目的とした。
**ライブラリ**: CyBase([https://cybase.org.au/](https://cybase.org.au/))、KNOTTIN([https://www.dsimb.inserm.fr/KNOTTIN/](https://www.dsimb.inserm.fr/KNOTTIN/))
**スクリーニングツール**:
HADDOCK 2.4
Kd, ΔG予測：PRODIGY webserver（[https://rascar.science.uu.nl/prodigy/](https://rascar.science.uu.nl/prodigy/)）
**MD**: GROMACS 2020.4を使用した200nsのMDシミュレーション
**評価法**:
- RMSD, RMSF, RoG, SASA
- 水素結合解析
- アラニンスキャン
- MM/PBSA, per-residue energy decomposition, PCA, FEL analysis

## 14. **Identification of novel PAD4 inhibitors using pharmacophore-based virtual screening, molecular docking, and molecular dynamics simulations**

https://chemrxiv.org/engage/chemrxiv/article-details/670bd91ecec5d6c1420546ab

**出版年**: 2024
**IF**: なし　ChemRxiv
**標的**: PAD4 (PDB: 8GOD)
**標的選定理由**: PAD4の過剰発現は関節リウマチなどの自己免疫疾患に関連しており、PAD4を抑制することにより、炎症を抑える治療法として有望であるため。
**ライブラリ**: PubChem、ZINC、MolPort
**スクリーニング**: SwissADME→PharmIT [https://pharmit.csb.pitt.edu/](https://pharmit.csb.pitt.edu/)
**分子ドッキング**: MOE 2024 (有償)
**MD**: Desmond Schrödinger
**パラメータ**：TIP3P, OPLA -AA, 300K, 1 atm
**評価法**: RMSD, RMSF, Rg, contact analysis, MM-GBSA
