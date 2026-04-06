---
title: "📰 RDF 解析"
free: false
---

Figure7のRDF plotを作成します。

タンパク質－リガンドの相互作用におけるRDF（Radial Distribution Function、放射分布関数）は、リガンドの原子とタンパク質の特定の原子や残基との間の相互作用を距離に基づいて解析する重要な指標です。RDFプロットを使うことで、リガンドとタンパク質の結合部位の相互作用の強さや距離依存性を視覚的に理解することができます。

上記のSer239とArg399に当たるアミノ酸とリガンドとのRDF解析をしていきます。

まずindexファイルを作成します。

```
gmx make_ndx -f md_0_10.gro -o index.ndx
```

以下のように出力されます。

```powershell
Command line:
  gmx make_ndx -f md_0_10.gro -o index.ndx

Reading structure file
Going to read 0 old index file(s)
Analysing residue names:
There are:   890    Protein residues
There are:     1      Other residues
There are: 51780      Water residues
There are:    23        Ion residues
Analysing Protein...
Analysing residues not classified as Protein/DNA/RNA/Water and splitting into groups...

  0 System              : 168557 atoms
  1 Protein             : 13165 atoms
  2 Protein-H           :  6560 atoms
  3 C-alpha             :   890 atoms
  4 Backbone            :  2670 atoms
  5 MainChain           :  3562 atoms
  6 MainChain+Cb        :  4342 atoms
  7 MainChain+H         :  4414 atoms
  8 SideChain           :  8751 atoms
  9 SideChain-H         :  2998 atoms
 10 Prot-Masses         : 13165 atoms
 11 non-Protein         : 155392 atoms
 12 Other               :    29 atoms
 13 UNL                 :    29 atoms
 14 NA                  :    23 atoms
 15 Water               : 155340 atoms
 16 SOL                 : 155340 atoms
 17 non-Water           : 13217 atoms
 18 Ion                 :    23 atoms
 19 Water_and_ions      : 155363 atoms

 nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index

>
```

Ser239は550番目のアミノ酸なので、

`ri 550`と記入し、Enter,

`name 20 Ser239`と記入し、Enterを押してください。

続いて、

Arg399は710番目のアミノ酸に当たるので、

`ri 710` と記入し、Enter

`name 21 Arg399`と記入し、Enterを押してください。

`q` で戻ってください。

続いて、以下のコードを実行してください。

```powershell
gmx rdf -s md_0_10.tpr -f md_0_10_center.xtc -n index.ndx -o rdf_Ser239_Arg399.xvg -b 0 -e 10000 -bin 0.01 -tu ns  -norm rdf
```

- **`gmx rdf`**: GROMACSでRDFを計算するコマンド。
- **`b 0`**: シミュレーションの最初のフレームを使って解析を開始します（0 nsから開始）。
- **`e 10000`**: 解析を終了するフレーム（ここでは10,000 ns）。
- **`bin 0.01`**: RDFを計算する際のビン幅。0.01 nmごとに距離を区切ってRDFを計算します。
- **`tu ns`**: 時間単位をナノ秒（ns）に設定。
- **`norm rdf`**: 正規化方法を標準的なRDF（全体の平均密度で正規化）に指定。

以下のように出力されます。

```powershell
                        :-) GROMACS - gmx rdf, 2024 (-:

Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /home/shizuku/protein_ligand_complex_ver6
Command line:
  gmx rdf -s md_0_10.tpr -f md_0_10_center.xtc -n index.ndx -o rdf_Ser239_Arg399.xvg -b 0 -e 10000 -bin 0.01 -tu ns -norm rdf

Available static index groups:
 Group  0 "System" (168557 atoms)
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
 Group 11 "non-Protein" (155392 atoms)
 Group 12 "Other" (29 atoms)
 Group 13 "UNL" (29 atoms)
 Group 14 "NA" (23 atoms)
 Group 15 "Water" (155340 atoms)
 Group 16 "SOL" (155340 atoms)
 Group 17 "non-Water" (13217 atoms)
 Group 18 "Ion" (23 atoms)
 Group 19 "Water_and_ions" (155363 atoms)
 Group 20 "Ser239" (17 atoms)
 Group 21 "Arg399" (15 atoms)
Specify a selection for option 'ref'
(Reference selection for RDF computation):
(one per line, <enter> for status/groups, 'help' for help)
>
```

13を選択してください。

続いて、以下が出力されます。

```powershell
Available static index groups:
 Group  0 "System" (168557 atoms)
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
 Group 11 "non-Protein" (155392 atoms)
 Group 12 "Other" (29 atoms)
 Group 13 "UNL" (29 atoms)
 Group 14 "NA" (23 atoms)
 Group 15 "Water" (155340 atoms)
 Group 16 "SOL" (155340 atoms)
 Group 17 "non-Water" (13217 atoms)
 Group 18 "Ion" (23 atoms)
 Group 19 "Water_and_ions" (155363 atoms)
 Group 20 "Ser239" (17 atoms)
 Group 21 "Arg399" (15 atoms)
Specify any number of selections for option 'sel'
(Selections to compute RDFs for from the reference):
(one per line, <enter> for status/groups, 'help' for help, Ctrl-D to en
```

`20`を選択し、Enter、

`21`を選択し、Enter

ctrl-Dで終了します。終了すると、計算が始まります。

`rdf_Ser239_Arg399.xvg` が生成しているので、グラフにしてみてください。

![](https://storage.googleapis.com/zenn-user-upload/0611dc63811c-20241125.png)


:::message
本記事ではリガンド全体を指定しているので、リガンド全ての原子のDistanceが加算されているようです。ですので、横軸Distanceが通常より多く出ています。
なので、本論文よりも多めに出てしまっていることをご了承ください。
論文と同様にリガンドの酸素原子のみをindexファイルで特異的に指定できますので、ご興味ある方はそちらでお試しください。
:::

 **青線（Ser(239)）のピーク**

- **ピーク位置**: 青線のピークは **約4.0 nm** の距離に位置しています。このピークはSer(239)に対する参照粒子が、約4.0 nmの距離で最も多く存在していることを示しています。
- **g(r)の値**: g(r)の最大値は **約8.0** です。これは、Ser(239)の周りに4.0 nmの距離で、全体の平均密度に対して8倍の密度で粒子が集まっていることを示しています。これは、Ser(239)と他の粒子が4.0 nmの距離で強い相互作用を持っている可能性を示唆します。

**赤線（Arg(399)）のピーク**

- **ピーク位置**: 赤線のピークは **約4.8 nm** に位置しています。Arg(399)との相互作用がこの距離で最も強く、粒子が集まっていることを示しています。
- **g(r)の値**: Arg(399)に対するg(r)の最大値は **約6.0** です。Ser(239)に比べると、少し密度が低いことを示していますが、それでも強い相互作用を持っていることがわかります。

**2つのピークの比較**

- **距離の違い**: Ser(239)は **約4.0 nm**、Arg(399)は **約4.8 nm** の距離で最も強い相互作用を示しています。このことから、これら2つの残基が異なる距離でリガンドや他の分子と相互作用していることが示されています。
    - これは、リガンドや他の分子が特定の距離でこれらの残基に結合する、もしくは相互作用することを示唆します。
- **強度の違い**: g(r)のピーク値から、Ser(239)の方がArg(399)よりもやや強い相互作用を示しています。これは、Ser(239)周囲の局所密度が高いことを意味しています。

以上でRMSD, RMSF, RoG, hydrogen bond, RDF解析は終了です！
