---
title: "📰 Indexファイルの作成"
free: false
---

それでは実際に解析していきます。

RMSDではタンパク質のリガンドの重原子の位置がどれくらい変化するかを計測します。

リガンドの重原子については、デフォルト設定では選択できないので、それを選択できるようにします。以下を実行してください。

```powershell
gmx make_ndx -f md_0_10.gro -n index.ndx
```

`gmx make_ndx` コマンドは、GROMACSでインデックスグループを作成または編集するために使います。以下はそのコマンドの使い方です。

- `f em.gro`: 構造ファイル（ここでは`em.gro`）を指定します。このファイルからインデックスグループを生成します。
- `n index.ndx`: 既存のインデックスファイルを読み込んで編集する場合に使います。ここでは `index.ndx` を読み込んで編集できます。

以下のように出力されます。

```powershell
                      :-) GROMACS - gmx make_ndx, 2024 (-:

Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /home/shizuku/protein_ligand_complex_ver3
Command line:
  gmx make_ndx -f em.gro -n index.ndx

Reading structure file
Going to read 1 old index file(s)

  0 System              : 169799 atoms
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
 11 non-Protein         : 156634 atoms
 12 Other               :    29 atoms
 13 UNL                 :    29 atoms
 14 NA                  :    23 atoms
 15 Water               : 156582 atoms
 16 SOL                 : 156582 atoms
 17 non-Water           : 13217 atoms
 18 Ion                 :    23 atoms
 19 Water_and_ions      : 156605 atoms
 20 Protein_UNL         : 13194 atoms

 nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index

> 
```

今回新しくリガンドの重原子を作成したいため、

`13 & a H*`　と打ち、Enterしてください。

次に以下を打ち、21番目にリガンドの重原子を選択できるようにします。

`name 21 UNL_Heaby`

最後に`q`を打ち、インタラクティブモードから出てください。
