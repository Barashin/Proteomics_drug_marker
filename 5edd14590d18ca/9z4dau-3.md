---
title: "📰 RMSF解析"
free: false
---

RMSF（Root Mean Square Fluctuation）は、分子動力学シミュレーションで各原子や残基の時間に伴う変動量を測定する方法です。具体的には、各原子がシミュレーション中に平均位置からどれくらい揺れ動いたかを示します。これにより、分子のどの部分が柔軟で、どの部分が安定しているかがわかります。

```powershell
gmx rmsf -s md_0_10.tpr -f md_0_10_center.xtc -o rmsf.xvg -res
```

- **`gmx rmsf`**: GROMACSのRMSF解析を行うコマンド。
- **`s topol.tpr`**: シミュレーションのトポロジーファイル（`.tpr`）を指定します。
- **`f center.xtc`**: トラジェクトリーファイル（`.xtc`）を指定します。このファイルにはシミュレーションの座標データが入っています。
- **`o rmsf.xvg`**: RMSFの結果を出力するファイル（`.xvg`形式）を指定します。
- **`res`**: 残基ごとにRMSFを計算するオプションです。

以下のように出力されます。

```powershell
(base) shizuku@DESKTOP-5I5GHRA:~/protein_ligand_complex_ver2$ gmx rmsf -s md_0_10.tpr -f md_0_10_center.xtc -o rmsf.xvg -res
                        :-) GROMACS - gmx rmsf, 2024 (-:

Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /home/shizuku/protein_ligand_complex_ver2
Command line:
  gmx rmsf -s md_0_10.tpr -f md_0_10_center.xtc -o rmsf.xvg -res

Reading file md_0_10.tpr, VERSION 2024 (single precision)
Reading file md_0_10.tpr, VERSION 2024 (single precision)
Select group(s) for root mean square calculation
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
Select a group: 
```

3を選択してください。

`rmsf.xvg`が出てくるので、これをエクセルで書いてください。

以下のようになります。

![](https://storage.googleapis.com/zenn-user-upload/cd0512c896f8-20241125.png)

このグラフは、タンパク質の残基ごとの**RMSF (Root Mean Square Fluctuation)** を示しています。RMSFは、各残基の平均的な位置からの振動（フレキシビリティ）を評価する指標です。以下のように解釈できます：

- **X軸 (Residue)**: 残基番号を表しており、312番目の残基から開始しています。グラフに表示されている範囲は、312番目から約700番目の残基に対応します。
- **Y軸 (RMSF)**: 各残基の振動の大きさ（Å単位）を示しています。数値が大きいほど、その残基がシミュレーション中に大きく動いている、すなわちフレキシブルであることを示します。
- グラフの左側（残基番号が小さい部分）では、RMSFが比較的低く（1-2Å程度）、この領域の残基はあまり動いていない（安定している）ことを示しています。
- 逆に、グラフの右側（残基番号が大きい部分）では、RMSFが高く（最大4Å以上）、この領域の残基は比較的動きやすいことを意味します。
- 全体的に、特定の残基が局所的に大きなフレキシビリティを示している部分（例えば600番目付近）が見られ、この領域がタンパク質の可動性が高い部分であることを示唆しています。

これにより、動きが少ない安定した領域と、動きが大きく柔軟性が高い領域を特定することができます。リガンドの結合したタンパク質とリガンドが結合していないタンパク質をRMSF解析で比べて、リガンドの結合による影響をみることが多いです。
