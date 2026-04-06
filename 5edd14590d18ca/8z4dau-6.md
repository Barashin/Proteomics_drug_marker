---
title: "📰 リガンド-タンパク質の設定"
free: false
---

続いて、以下を実行してください。

```powershell
gmx make_ndx -f ligand_fix.acpype/ligand_fix_GMX.gro -o index_ligand_fix_GMX.ndx
```

このコマンドを使うと、特定の原子や分子のグループ（インデックス）を定義し、後の解析やシミュレーションで使用することができます。

- `gmx make_ndx`: GROMACSのインデックス作成・編集ツールを起動します。
- `f ligand_fix_GMX.gro`: インデックス作成の対象となる構造ファイル（この場合、`ligand_fix_GMX.gro`）を指定します。この構造ファイルに基づいて、特定の原子や分子をインデックスとして定義します。
- `o index_ligand_fix_GMX.ndx`: 出力するインデックスファイルの名前を指定します。結果として作成されるインデックスファイルは`index_ligand_fix_GMX.ndx`という名前になります。

以下のように出力されます。

```powershell
                      :-) GROMACS - gmx make_ndx, 2024 (-:

Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /home/shizuku/protein_ligand_complex_ver2
Command line:
  gmx make_ndx -f ligand_GMX.gro -o index_ligand_GMX.ndx

Reading structure file
Going to read 0 old index file(s)
Analysing residue names:
There are:     1      Other residues
Analysing residues not classified as Protein/DNA/RNA/Water and splitting into groups...

  0 System              :    29 atoms
  1 Other               :    29 atoms
  2 UNL                 :    29 atoms

 nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index
```

 `0 & ! a H*`　を入力してください。

- `0 & ! a H*` は、「システム全体から水素原子を除いたグループ」を作成するコマンドです。
- 水素原子を除外したインデックスを使って、特定の解析やシミュレーションで水素以外の原子にフォーカスできます。

以下のように出力されます。

```powershell
Copied index group 0 'System'
Found 14 atoms with name H*
Complemented group: 15 atoms
Merged two groups with AND: 29 15 -> 15
```

`q`　と打ち、対話から出てください。

続いて、以下を実行してください。

```powershell
gmx genrestr -f ligand_fix.acpype/ligand_fix_GMX.gro -n index_ligand_fix_GMX.ndx -o posre_ligand_fix_GMX.itp -fc 1000 1000 1000
```

このコマンドは、GROMACSで位置拘束（Position Restraint）を作成するために使用されます。位置拘束は、シミュレーション中に特定の原子や分子の位置をある程度固定したい場合に便利です。

- **`gmx genrestr`**:
    - GROMACSの位置拘束ファイル（`posre.itp` ファイル）を生成するためのコマンドです。このファイルは、指定された原子や分子の位置を固定（または拘束）し、シミュレーション中の特定の動きを制限するために使用されます。
- **`f ligand_fix_GMX.gro`**:
    - `ligand_fix_GMX.gro` という構造ファイルを入力ファイルとして使用します。このファイル内の原子座標が位置拘束の対象となります。
- **`n index_ligand_fix_GMX.ndx`**:
    - `index_ligand_fix_GMX.ndx` というインデックスファイルを指定します。このファイルは、拘束をかけたい特定の原子グループを定義します。例えば、特定のリガンドや特定の原子のみを拘束したい場合に、このインデックスファイルを使用します。
- **`o posre_ligand_GMX.itp`**:
    - `posre_ligand_GMX.itp` という出力ファイル名を指定します。生成されるファイルは、位置拘束のパラメータが記述された `.itp` ファイルです。このファイルは、後でシミュレーションで使用するトポロジーファイルにインクルードされます。
- **`fc 1000 1000 1000`**:
    - `fc` は、**拘束力定数（force constants）**を指定します。この場合、x軸、y軸、z軸の各方向において1000 kJ/mol/nm²の拘束力定数を設定します。この値が大きいほど、原子の位置を強く拘束します。
    - `1000 1000 1000` は、各軸ごとの拘束力定数を表しており、この例ではすべての方向（x, y, z）で同じ拘束力がかかるように設定されています。

以下のように出力されます。

```powershell
(base) shizuku@DESKTOP-5I5GHRA:~/protein_ligand_complex_ver3$ gmx genrestr -f ligand_fix_GMX.gro -n index_ligand_fix_GMX.ndx -o posre_ligand_GMX.itp -fc 1000 1000 1000
                      :-) GROMACS - gmx genrestr, 2024 (-:

Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /home/shizuku/protein_ligand_complex_ver3
Command line:
  gmx genrestr -f ligand_fix_GMX.gro -n index_ligand_fix_GMX.ndx -o posre_ligand_GMX.itp -fc 1000 1000 1000

Reading structure file
Select group to position restrain
Group     0 (         System) has    29 elements
Group     1 (          Other) has    29 elements
Group     2 (            UNL) has    29 elements
Group     3 (   System_&_!H*) has    15 elements
Select a group:
```

3を選択してください。

続いて、Topol.topに以下の**Ligand position restraints**を追加してください。

```powershell
;
;	File 'topol.top' was generated
;	By user: unknown (1000)
;	On host: DESKTOP-5I5GHRA
;	At date: Thu Oct  3 20:12:32 2024
;
;	This is a standalone topology file
;
;	Created by:
;	                     :-) GROMACS - gmx pdb2gmx, 2024 (-:
;
;	Executable:   /usr/local/gromacs/bin/gmx
;	Data prefix:  /usr/local/gromacs
;	Working dir:  /home/shizuku/protein_ligand_complex_ver3
;	Command line:
;	  gmx pdb2gmx -f protein.pdb -o protein_processed.gro -ignh
;	Force field was read from the standard GROMACS share directory.
;
; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"

[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
 ca       ca          0.00000  0.00000   A     3.31521e-01   4.13379e-01 ; 1.86  0.0988
 ce       ce          0.00000  0.00000   A     3.31521e-01   4.13379e-01 ; 1.86  0.0988
 c2       c2          0.00000  0.00000   A     3.31521e-01   4.13379e-01 ; 1.86  0.0988
 c3       c3          0.00000  0.00000   A     3.39771e-01   4.51035e-01 ; 1.91  0.1078
 os       os          0.00000  0.00000   A     3.15610e-01   3.03758e-01 ; 1.77  0.0726
 oh       oh          0.00000  0.00000   A     3.24287e-01   3.89112e-01 ; 1.82  0.0930
 ho       ho          0.00000  0.00000   A     5.37925e-02   1.96648e-02 ; 0.30  0.0047
 ha       ha          0.00000  0.00000   A     2.62548e-01   6.73624e-02 ; 1.47  0.0161
 h1       h1          0.00000  0.00000   A     2.42200e-01   8.70272e-02 ; 1.36  0.0208
 hc       hc          0.00000  0.00000   A     2.60018e-01   8.70272e-02 ; 1.46  0.0208

; Include chain topologies
#include "topol_Protein_chain_A.itp"
#include "topol_Protein_chain_B.itp"

; Include ligand topology
#include "ligand_fix.acpype/ligand_fix_GMX.itp"

; Ligand position restraints (←ここを追加)
#ifdef POSRES
#include "posre_ligand_fix_GMX.itp"
#endif

; Include water topology
#include "amber99sb-ildn.ff/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "amber99sb-ildn.ff/ions.itp"

[ system ]
; Name
3-PHOSPHOSHIKIMATE 1-CARBOXYVINYLTRANSFERASE

[ molecules ]
; Compound        #mols
Protein_chain_A     1
Protein_chain_B     1
UNL		    1
SOL         52194
NA               23

```

GROMACSのシミュレーションにおいて、特定の分子（タンパク質、リガンドなど）を動きにくくするために位置拘束を使います。これは、特に平衡化（equilibration）の段階で、特定の部分を安定化させたい場合に有効です。

上記の設定で、タンパク質とリガンドの両方を**同時に拘束**しています。リガンドとタンパク質を別々に拘束したいときは[GROMACSのチュートリアル](http://www.mdtutorials.com/gmx/complex/06_equil.html)を参考にしてください。

続いて以下を実行してください。

```powershell
gmx make_ndx -f em.gro -o index.ndx
```

`gmx make_ndx -f em.gro -o index.ndx` は、GROMACSでインデックスファイル（`index.ndx`）を作成・編集するためのコマンドです。このコマンドを使うことで、特定の分子や原子をグループ化し、そのグループをシミュレーションの解析や操作で使用できるようになります。

```powershell
Reading structure file
Going to read 0 old index file(s)
Analysing residue names:
There are:   890    Protein residues
There are:     1      Other residues
There are: 52194      Water residues
There are:    23        Ion residues
Analysing Protein...
Analysing residues not classified as Protein/DNA/RNA/Water and splitting into groups...

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

 nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index
```

 `1 | 13`　と書き、Enterを押してください。

以下のように出力されます。

```powershell
Copied index group 1 'Protein'
Copied index group 13 'UNL'
Merged two groups with OR: 13165 29 -> 13194
```

`q`　と書いて、実行してください。

この操作により、**タンパク質とリガンドが1つのグループとして扱えるようになります**。
