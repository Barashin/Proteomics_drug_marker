---
title: "📰 RMSD解析"
free: false
---

RMSD（Root Mean Square Deviation）は、分子動力学シミュレーションにおいて、分子の構造変化を時間に渡って追跡するための解析手法です。RMSDは、シミュレーション中の各フレームの原子の位置を、参照構造と比較して、どれだけ変化しているかを測定します。

今回は先ほど作成した `UNL_Heaby`　がシミュレーション中でどのくらい変化しているかを見ます。

まずはcdでMD simulationで得られたデータが入っているディレクトリに移動してください。

`md_0_10.tpr`、　`md_0_10_center.xtc` 、`index.ndx`の三つが入っているか解説してください。

![](https://storage.googleapis.com/zenn-user-upload/2c404e29f4cf-20241125.jpg)

RMSD解析

```powershell
gmx rms -s md_0_10.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd.xvg
```

このコマンドは、GROMACSの`gmx rms`ツールを使って、分子動力学（MD）シミュレーションのRMSD（Root Mean Square Deviation、二乗平均平方根偏差）を計算するためのものです。特定の構造（例: エネルギー最小化後の構造）とシミュレーション中の構造との違いを時間経過に伴って測定します。RMSDは、シミュレーションの安定性や系の構造変化を評価するための重要な指標です。

gmx rmsは、GROMACSのRMSDを計算するためのコマンドです。このツールは、リファレンス構造（通常はエネルギー最小化後の構造など）に対する各時点のシミュレーション構造の座標を比較し、RMSDを計算します。

- `gmx rms`は、GROMACSのRMSDを計算するためのコマンドです。このツールは、リファレンス構造（通常はエネルギー最小化後の構造など）に対する各時点のシミュレーション構造の座標を比較し、RMSDを計算します。
- `s`はリファレンス構造を指定するオプションです。ここでは、エネルギー最小化後の構造情報が入った`em.tpr`ファイルを指定しています。このファイルは、RMSD計算時の基準となる構造です。通常はシミュレーション開始時の最適化された構造を使用します。
- `f`は、RMSDを計算するために使用するシミュレーションの軌道（トラジェクトリ）ファイルを指定します。ここでは、軌道ファイル`md_0_10_center.xtc`を指定しています。このファイルにはシミュレーション中の構造の座標情報が含まれており、特定の基準構造との比較に使われます。
- `n`はインデックスファイルを指定するオプションです。`index.ndx`ファイルには、特定の原子や分子グループが定義されています。このファイルを使うことで、RMSDを計算する際にどの部分（タンパク質全体、特定のサブユニット、リガンドなど）を比較対象とするかを指定できます。たとえば、タンパク質の主鎖のみや、リガンドのRMSDを計算したい場合に便利です。
- `tu`は時間単位を指定するオプションで、ここでは`ns`（ナノ秒）を指定しています。これにより、出力されるRMSDデータの時間単位がナノ秒になります。デフォルトではピコ秒（ps）ですが、ナノ秒の方が長時間シミュレーションでは扱いやすいです。
- `o`は、RMSD計算結果の出力ファイル名を指定します。ここでは、`rmsd.xvg`という名前のファイルにRMSDのデータが出力されます。このファイルは、Xmgraceなどのプロットツールで可視化できる形式（`.xvg`）です。時間とRMSD値の変化をグラフ化する際に使用します。

以下のように出力されます。

```powershell
Command line:
  gmx rms -s em.tpr -f md_0_10_center.xtc -n index.ndx -tu ns -o rmsd.xvg

Reading file em.tpr, VERSION 2024 (single precision)
Reading file em.tpr, VERSION 2024 (single precision)
Select group for least squares fit
Group     0 (         System) has 169799 elements
Group     1 (        Protein) has 13165 elements
Group     2 (      Protein-H) has  6560 elements
Group     3 (        C-alpha) has   890 elements
Group     4 (       Backbone) has  2670 elements
Group     5 (      MainChain) has  3562 elements
Group     6 (   MainChain+Cb) has  4342 elements
Group     7 (    MainChain+H) has  4414 elements
Group     8 (      SideChain) has  8751 elements
Group     9 (    SideChain-H) has  2998 elements
Group    10 (    Prot-Masses) has 13165 elements
Group    11 (    non-Protein) has 156634 elements
Group    12 (          Other) has    29 elements
Group    13 (            UNL) has    29 elements
Group    14 (             NA) has    23 elements
Group    15 (          Water) has 156582 elements
Group    16 (            SOL) has 156582 elements
Group    17 (      non-Water) has 13217 elements
Group    18 (            Ion) has    23 elements
Group    19 ( Water_and_ions) has 156605 elements
Group    20 (    Protein_UNL) has 13194 elements
Group    21 (      UNL_Heaby) has    14 elements
Select a group: 
```

4を選択してください。

```powershell
Select group for RMSD calculation
Group     0 (         System) has 169799 elements
Group     1 (        Protein) has 13165 elements
Group     2 (      Protein-H) has  6560 elements
Group     3 (        C-alpha) has   890 elements
Group     4 (       Backbone) has  2670 elements
Group     5 (      MainChain) has  3562 elements
Group     6 (   MainChain+Cb) has  4342 elements
Group     7 (    MainChain+H) has  4414 elements
Group     8 (      SideChain) has  8751 elements
Group     9 (    SideChain-H) has  2998 elements
Group    10 (    Prot-Masses) has 13165 elements
Group    11 (    non-Protein) has 156634 elements
Group    12 (          Other) has    29 elements
Group    13 (            UNL) has    29 elements
Group    14 (             NA) has    23 elements
Group    15 (          Water) has 156582 elements
Group    16 (            SOL) has 156582 elements
Group    17 (      non-Water) has 13217 elements
Group    18 (            Ion) has    23 elements
Group    19 ( Water_and_ions) has 156605 elements
Group    20 (    Protein_UNL) has 13194 elements
Group    21 (      UNL_Heaby) has    14 elements
Select a group: 
```

21を選択してください。

上記の設定は、シミュレーション中にリガンド（UNL_Heavy）が**タンパク質の主鎖構造に対してどのように動いたか**、その相対的な変化を計算するためのものです。RMSD値が時間とともに安定している場合、リガンドは基準構造（タンパク質の主鎖）に対して安定していることを意味し、大きく変動する場合はリガンドがシミュレーション中に動き回ったことを示します。

`ls`コマンドを押して、以下の`rmsd.xvg`あるか確認してください。

`rmsd.xvg`をテキストで開いた後、エクセルなどで張り付けてグラフを作成してください。

以下のようにRMSDができます。

![](https://storage.googleapis.com/zenn-user-upload/6d1263c74187-20241125.png)

このRMSDの変動パターンから、シミュレーションの最初の25 nsでは大きな構造変化が起こっていたが、以降は徐々に安定し、最終的には構造が収束していると考えられます。これにより、リガンドが結合部位に安定した可能性が示唆されますが、いくつかの時間点で再び変動が見られるため、完全に安定した結合ではない可能性もあります。
