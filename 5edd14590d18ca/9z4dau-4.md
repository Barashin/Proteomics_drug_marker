---
title: "📰 RoG解析"
free: false
---

**RoG解析**（Radius of Gyration、回転半径）は、分子がシミュレーション中にどれだけコンパクトにまとまっているかを測定する指標です。これは、分子全体の重心から各原子の距離を考慮して計算され、タンパク質や分子の折りたたみ状態や構造の安定性を評価するのに使われます。

```powershell
gmx gyrate -s md_0_10.tpr -f md_0_10_center.xtc -o gyrate.xvg
```

コマンド `gmx gyrate` は、分子の重心からの距離を基準に、シミュレーション中のタンパク質や分子の回転半径（gyration radius）を計算するために使用されます。

以下のように出力されます。

```powershell
Command line:
  gmx gyrate -s md_0_10.tpr -f md_0_10_center.xtc -o gyrate.xvg

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
Specify a selection for option 'sel'
(Select group to compute gyrate radius):
(one per line, <enter> for status/groups, 'help' for help)
>
```

1を選択してください。

`gyrate.xvg`　が出力されます。一行目と二行目のデータを使って、散布図を書いてください。

以下のように書けます。

![](https://storage.googleapis.com/zenn-user-upload/b4178ba6b413-20241125.png)

- シミュレーションの最初の部分では、RoGが35Å付近で変動していますが、徐々に減少し、約30Åで安定しています。
- 250 ns以降は、RoGが一定の範囲で安定しており、これはシミュレーションが収束し、タンパク質の全体的な構造が安定したことを示しています。
