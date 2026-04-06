---
title: "📰 後処理"
free: false
---

周期的な境界条件下で行われるシミュレーションと同様に、分子は "壊れて "表示されたり、ボックスを "飛び跳ね "て行ったり来たりすることがあります。 これを解決するために、タンパク質を再中心化します。

以下のコードを実行してください。

```powershell
gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
```

- `gmx trjconv`: これはGROMACSのツールで、軌跡ファイルを変換したり処理するために使われる `trjconv` コマンドです。
- `s md_0_10.tpr`: これは構造情報が含まれた入力ファイル（`.tpr` ファイル）です。このファイルは、軌跡データの処理に必要なトポロジーや他の情報を提供します。
- `f md_0_10.xtc`: これはシミュレーションの軌跡ファイル（`.xtc` ファイル）です。シミュレーション中の原子の座標が含まれています。
- `o md_0_10_center.xtc`: これは出力ファイルの指定です。変換された軌跡を保存するための新しいファイルで、この場合は `md_0_10_center.xtc` として出力されます。
- `center`: 指定した分子やグループの重心をシミュレーションボックスの中心に移動させます。これにより、特定の分子が常にボックスの中心付近に維持されるようになります。
- `pbc mol`: `pbc` は周期境界条件（PBC: Periodic Boundary Conditions）を扱うオプションです。この場合、`mol` は分子単位で周期境界条件を適用し、分子がボックスの外に出た場合も分子全体を元の位置に戻します（分子が切り裂かれるのを防ぎます）。
- `ur compact`: `ur` は "unit cell representation" を指定します。このオプションで `compact` を指定すると、最もコンパクトな形で分子を配置するようにします。これは視覚化や解析の際に便利です。

以下のように出力されます。

```powershell
Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /home/shizuku/protein_ligand_complex_ver3
Command line:
  gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact

Note that major changes are planned in future for trjconv, to improve usability and utility.
Will write xtc: Compressed trajectory (portable xdr format): xtc
Reading file md_0_10.tpr, VERSION 2024 (single precision)
Reading file md_0_10.tpr, VERSION 2024 (single precision)
Select group for centering
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

1を選択してください

```powershell
Select group for output
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

０を選択してください。

これで後処理は終了です！

**※これでもうまくいかない場合は、以下を実行してください。**
簡易的にコードを示しておきます。
詳しくはAg_smith氏のQiita記事「[失敗しないGROMACSトラジェクトリの変換方法](https://qiita.com/Ag_smith/items/d1984cd11bb24919e5d5)」を参考してください。


```
gmx make_ndx -f md_0_10.gro -n index.ndx
```

上記でProtein-ligandのcomplexの選択肢を作ります。

```
gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_whole.xtc -n index.ndx -pbc whole
```

0を選択します。

```
gmx trjconv -s md_0_10.tpr -f md_0_10_whole.xtc -o md_0_10_whole_cluster.xtc -n index.ndx -pbc cluster
```

21を選択します。

0を選択します。

```
gmx trjconv -s md_0_10.tpr -f md_0_10_whole_cluster.xtc -o md_0_10_whole_cluster_center.xtc -n index.ndx -pbc mol -ur compact -center
```

21を選択します。

0を選択します。

```
gmx trjconv -s md_0_10.tpr -f md_0_10_whole_cluster_center.xtc -o md_0_10_complete.xtc -n index.ndx -center -fit rot+trans
```

21を選択します。

21を選択します。

0を選択します。


`md_0_10_complete.xtc`をファイルとして、PyMOLで見るとうまくMDができていると思います！
