---
title: "📰 水素結合数測定"
free: false
---

MDシミュレーションにおける水素結合数の測定は、分子間および分子内での結合の強さや安定性を評価する重要な解析手法です。水素結合（H-bond）は、分子が互いに相互作用する際に形成される強力な非共有結合の一種で、タンパク質-リガンド相互作用や分子の安定性において非常に重要です。

以下を実行して、タンパク質とリガンドの水素結合の数を計算してください。

※`gmx hbond` コマンドでは[うまく解析できないことがあるそうです](https://gromacs.bioexcel.eu/t/inconsistency-of-gmx-hbond-analysis/8785)。

```powershell
gmx hbond-legacy -s md_0_10.tpr -f md_0_10_center.xtc -num hydrogen-bonds.xvg -tu ns
```

`gmx hbond-legacy` コマンドを使用して、水素結合を解析し、`-num` オプションで結合の数を `hydrogen-bonds.xvg` に出力するコマンドの説明に移ります。

以下のように出力されます。

```powershell
Command line:
  gmx hbond -s md_0_10.tpr -f md_0_10_center.xtc -num hydrogen-bonds.xvg -tu ns

Reading file md_0_10.tpr, VERSION 2024 (single precision)
Reading file md_0_10.tpr, VERSION 2024 (single precision)
Available static index groups:
 Group  0 "System" (169799 atoms)
 Group  1 "Protein" (13165 atoms)
 Group  2 "Protein-H" (6560 atoms)
 Group  3 "C-alpha" (890 atoms)
 Group  4 "Backbone" (2670 atoms)
 Group  5 "MainChain" (3562 atoms)
 Group  6 "MainChain+Cb" (4342 atoms)
 Group  7 "MainChain+H" (4414 atoms)
 Group  8 "SideChain" (8751 atoms)
 Group  9 "SideChain-H" (2998 atoms)
 Group 10 "Prot-Masses" (13165 atoms)
 Group 11 "non-Protein" (156634 atoms)
 Group 12 "Other" (29 atoms)
 Group 13 "UNL" (29 atoms)
 Group 14 "NA" (23 atoms)
 Group 15 "Water" (156582 atoms)
 Group 16 "SOL" (156582 atoms)
 Group 17 "non-Water" (13217 atoms)
 Group 18 "Ion" (23 atoms)
 Group 19 "Water_and_ions" (156605 atoms)
Specify a selection for option 'r'
(Reference selection, relative to which the search for hydrogen bonds in target selection will develop.):
(one per line, <enter> for status/groups, 'help' for help)
> 
```

1を選択してください。

```powershell
Available static index groups:
 Group  0 "System" (169799 atoms)
 Group  1 "Protein" (13165 atoms)
 Group  2 "Protein-H" (6560 atoms)
 Group  3 "C-alpha" (890 atoms)
 Group  4 "Backbone" (2670 atoms)
 Group  5 "MainChain" (3562 atoms)
 Group  6 "MainChain+Cb" (4342 atoms)
 Group  7 "MainChain+H" (4414 atoms)
 Group  8 "SideChain" (8751 atoms)
 Group  9 "SideChain-H" (2998 atoms)
 Group 10 "Prot-Masses" (13165 atoms)
 Group 11 "non-Protein" (156634 atoms)
 Group 12 "Other" (29 atoms)
 Group 13 "UNL" (29 atoms)
 Group 14 "NA" (23 atoms)
 Group 15 "Water" (156582 atoms)
 Group 16 "SOL" (156582 atoms)
 Group 17 "non-Water" (13217 atoms)
 Group 18 "Ion" (23 atoms)
 Group 19 "Water_and_ions" (156605 atoms)
Specify a selection for option 't'
(Target selection, relative to which the search for hydrogen bonds in reference selection will develop.):
(one per line, <enter> for status/groups, 'help' for help)
> 
```

13を選択してください。

`hydrogen-bonds.xvg` が出力されているので、こちらの一行目と二行目のデータを使ってグラフを書いてください。

![](https://storage.googleapis.com/zenn-user-upload/a5635ea0657f-20241125.png)

このグラフから、最初の25 nsではタンパク質とリガンドの間に比較的安定した水素結合が存在していますが、シミュレーションが進行するにつれて、水素結合の数が減少し、不安定な状態になっていることがわかります。これは、リガンドが結合ポケット内で再配置されているか、結合から外れ始めていることを示している可能性があります。
