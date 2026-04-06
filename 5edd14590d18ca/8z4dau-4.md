---
title: "📰 水とイオン追加"
free: false
---

システムを水溶液中に配置し、必要に応じて中和します。システムに12 Åのバッファ（水ボックスの大きさ）を持つ水分子を追加します。

1. **水ボックスの追加**：

    ```powershell
    gmx editconf -f complex.gro -o boxed.gro -c -d 1.2 -bt cubic
    ```

    このコマンドは12 Åのバッファを設定し、システムを水分子で囲みます。

2. **溶媒（水）の追加**：

    ```bash
    gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
    ```

3. **イオンの追加**（中和）：

ions.mdpを作成します。これについては、論文内で記述がないので、[gromacsチュートリアル](http://www.mdtutorials.com/gmx/complex/04_ions.html)から取ってきたものをそのまま使います。

```
; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
rlist		    = 1.0		; Cut-off for making neighbor list (short range forces)
coulombtype	    = cutoff	; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; long range electrostatic cut-off
rvdw		    = 1.0		; long range Van der Waals cut-off
pbc             = xyz 		; Periodic Boundary Conditions
```

以下を実行してください。

```powershell
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
```

以下のようなエラーが出るかもしれません。

```powershell
Fatal error:
Syntax error - File ligand_fix_GMX.itp, line 3
Last line read:
'[ atomtypes ]'
Invalid order for directive atomtypes

```

`ligand_fix_GMX.itp`ファイルの中身が正しい順序になっていないからだそうです。

[このリンク](https://gromacs.bioexcel.eu/t/invalid-order-for-directive-atomtypes-error/3859/2)を参考にすると、topol topに、ligand _fix_GMX.itpの[ atomtypes ]の部分を取ってきて、張り付けてください。

以下が修正したtopol.topです。

```
; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"

[ atomtypes ]（←ここ追加）
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
UNL		    　　　　　 1
SOL             52217

```

現在のligand_fix_GMX.itpは以下です。

一部は消して、[ moleculetype ]のところはUNLとしています。

（topol.topの[ molecules ]のところと名前を合わせましょう。）

```
; ligand_fix_GMX.itp created by acpype (v: 2023.11.14) on Thu Oct  3 21:04:33 2024

（ここから）
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
 （ここまでは消去）

[ moleculetype ]
;name            nrexcl
 UNL               3　 （←ここは変更）

[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
     1   ca     1   UNL     C    1    -0.148000     12.01000 ; qtot -0.148
     2   ca     1   UNL    C1    2    -0.102000     12.01000 ; qtot -0.250
     3   ca     1   UNL    C2    3    -0.147000     12.01000 ; qtot -0.397
     4   ca     1   UNL    C3    4    -0.011800     12.01000 ; qtot -0.409
     5   ca     1   UNL    C4    5    -0.159300     12.01000 ; qtot -0.568
     6   ca     1   UNL    C5    6     0.127100     12.01000 ; qtot -0.441
     7   ce     1   UNL    C6    7    -0.108200     12.01000 ; qtot -0.549
     8   c2     1   UNL    C7    8    -0.176200     12.01000 ; qtot -0.725
     9   c3     1   UNL    C8    9     0.142300     12.01000 ; qtot -0.583
    10   os     1   UNL     O   10    -0.411600     16.00000 ; qtot -0.995
    11   c3     1   UNL    C9   11     0.177700     12.01000 ; qtot -0.817
    12   c3     1   UNL   C10   12     0.151100     12.01000 ; qtot -0.666
    13   c3     1   UNL   C11   13    -0.089100     12.01000 ; qtot -0.755
    14   oh     1   UNL    O1   14    -0.595799     16.00000 ; qtot -1.351
    15   ho     1   UNL     H   15     0.399000      1.00800 ; qtot -0.952
    16   oh     1   UNL    O2   16    -0.508100     16.00000 ; qtot -1.460
    17   ho     1   UNL    H1   17     0.422000      1.00800 ; qtot -1.038
    18   ha     1   UNL    H2   18     0.151000      1.00800 ; qtot -0.887
    19   ha     1   UNL    H3   19     0.136000      1.00800 ; qtot -0.751
    20   ha     1   UNL    H4   20     0.135000      1.00800 ; qtot -0.616
    21   ha     1   UNL    H5   21     0.126000      1.00800 ; qtot -0.490
    22   ha     1   UNL    H6   22     0.139000      1.00800 ; qtot -0.351
    23   h1     1   UNL    H7   23     0.036700      1.00800 ; qtot -0.314
    24   h1     1   UNL    H8   24     0.055700      1.00800 ; qtot -0.258
    25   h1     1   UNL    H9   25     0.055700      1.00800 ; qtot -0.203
    26   h1     1   UNL   H10   26     0.049700      1.00800 ; qtot -0.153
    27   hc     1   UNL   H11   27     0.051033      1.00800 ; qtot -0.102
    28   hc     1   UNL   H12   28     0.051033      1.00800 ; qtot -0.051
    29   hc     1   UNL   H13   29     0.051033      1.00800 ; qtot -0.000

```

すると以下が実行できるようになります。

```powershell
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
```

このコマンドは、GROMACSのシミュレーション準備ツールである`grompp`を使用して、MDシミュレーションの実行準備を行います。具体的には、指定されたパラメータファイル（`ions.mdp`）、座標ファイル（`solvated.gro`）、トポロジーファイル（`topol.top`）を基に、実行可能なシミュレーションの入力ファイル（`ions.tpr`）を生成します。

続いて、以下を実行してください。

```powershell
gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral
```

以下のように表示されます。

```powershell
Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /home/shizuku/protein_ligand_complex_ver3
Command line:
  gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral

Reading file ions.tpr, VERSION 2024 (single precision)
Reading file ions.tpr, VERSION 2024 (single precision)
Will try to add 23 NA ions and 0 CL ions.
Select a continuous group of solvent molecules
Group     0 (         System) has 169845 elements
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
Group    11 (    non-Protein) has 156680 elements
Group    12 (          Other) has    29 elements
Group    13 (            UNL) has    29 elements
Group    14 (          Water) has 156651 elements
Group    15 (            SOL) has 156651 elements
Group    16 (      non-Water) has 13194 elements
Select a group:
```

15を選択してください。

`solvated_ions.gro`ファイルが生成されていると思います

これにて、水の追加と中和は終了です。
