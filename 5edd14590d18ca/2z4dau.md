---
title: "構成（目次）"
free: false
---

**まえがき**

**構成（目次）**

**📘 第0章：in silico創薬の流れと同様の流れで作成されている論文たち**

- 📰 本書のIn silico創薬論文の流れと概要
- 📰 同様の流れで作成されている論文たち

**📘 第1章：ターゲット準備、minimization**

- 📰 ターゲットの選定
- 📰 タンパク質の準備
- 📰 タンパク質の安定性評価

**📘 第2章：Homology modelingによる構造予測**

- 📰 Homology ModelingとSWISS-MODELとは
- 📰 SWISS-MODELを使ったHomology modeling
- 📰 Homology modelingの結果

**📘 第3章：化合物ライブラリの構築**

- 📰 CMNPDからの抽出
- 📰 FAFDrug4 serverによるフィルタリング

**📘 第4章：化合物ライブラリまとめ**

- 📰 本記事で紹介する化合物データベースの比較
- 📰 ZINC15
- 📰 PubChem
- 📰 ChEMBL
- 📰 COCONUT（COlleCtion of Open Natural prodUcTs）
- 📰 マイナーだが、よく論文で使用されているDB
- 📰 他ライブラリの検索

**📘 第5章：RDKitを使った化合物のフィルタリング**

- 📰 環境構築
- 📰 RDKitを用いた化合物のフィルタリング
- 📰 出力結果

**📘 第6章：In silico screening**

- 📰 PyRxによる化合物ライブラリのエネルギー最小化
- 📰 In silicoスクリーニング
- 📰 結果

**📘 第7章：Sminaを使ったin silicoスクリーニング**

- 📰 Sminaの説明とタンパク質準備
- 📰 化合物ライブラリのイオン原子の除去
- 📰 化合物ライブラリの水素、電荷付加、3D構造化、エネルギー最小化
- 📰 Grid Boxの設定
- 📰 Sminaの環境構築
- 📰 configファイルの設定
- 📰 Sminaによるin silico screening
- 📰 結果

**📘 第8章：タンパク質-リガンドの相互作用可視化**

- 📰 PLIPとは
- 📰 PLIPの環境構築
- 📰 Cathepsin Kと阻害剤の相互作用解析
- 📰 スクリーニングした後の化合物保存の際の注意点
- 📰 PoseViewとは
- 📰 PoseViewを使ったタンパク質ーリガンドの相互作用解析

**📘 第9章：MD simulationのための環境構築**

- 📰 WindowsへのLinuxの導入
- 📰 GROMACSのインストール
- 📰 MacでのGROMACSのインストール（検証中)

**📘 第10章：MD simulationを用いたタンパク質の安定化**

- 📰 タンパク質準備
- 📰 PDBからGROMACSトポロジーへの変換
- 📰 ボックス作成と溶媒の追加
- 📰 イオンの追加（中性化）
- 📰 エネルギー最小化
- 📰 平衡化（NVT & NPT)
- 📰 MDシミュレーション
- 📰 可視化
- 📰 RMSD解析
- 📰 安定な構造の抽出

**📘 第11章：リガンド-タンパク質 MDシミュレーション（分子動力学シミュレーション）**

- 📰 リガンドの準備
- 📰 タンパク質の準備
- 📰 リガンド-タンパク質の準備
- 📰 水とイオン追加
- 📰 エネルギー最小化
- 📰 リガンド-タンパク質の設定
- 📰 NVT（定温）とNPT（定圧定温）の平衡化
- 📰 MDシミュレーション
- 📰 後処理
- 📰 可視化

**📘 第12章：RMSD, RMSF, RoG, 水素結合数測定, RDF解析**

- 📰 Indexファイルの作成
- 📰 RMSD解析
- 📰 RMSF解析
- 📰 RoG解析
- 📰 水素結合数測定
- 📰 RDF 解析

**📘 第13章：MM-GBSA/PBSA、decomposition、アラニンスキャン解析**

- 📰 gmx_MMPBSAの環境構築
- 📰 MM-GBSA/PBSA、decomposition、アラニンスキャン解析の実行
- 📰 MM-GBSA/PBSA解析結果
- 📰 decomposition解析結果
- 📰 アラニンスキャン解析結果

**📘 第14章：Water swapping解析**

- 📰 Waterswapping解析とは
- 📰 OpenBioSimの環境構築
- 📰 GROMACSファイルをAMBERファイルに変換
- 📰 Waterswapping解析実行
- 📰 結果

**📘 第15章：物性評価**

- 📰 物性評価とは
- 📰 Swiss ADMEによる物性評価
- 📰 Deep-PKによる物性評価

**あとがき**

**In silico創薬レクチャー、コンサルティング、共同研究について**
