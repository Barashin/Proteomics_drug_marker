---
title: "📰 リガンド-タンパク質の準備"
free: false
---

リガンドとタンパク質を結合した複合体の`gro`ファイルを生成します。これには手動で`protein_processed.gro`と`ligand_fix_GMX.gro`を結合するか、以下の方法で新しい構造ファイルを作成します。

1. タンパク質とリガンドを手動で`complex.gro`ファイルにまとめます。

`protein_processed.gro`　に`ligand_fix_GMX.gro`　の値を加え、以下の部分を変更してください。

```
3-PHOSPHOSHIKIMATE 1-CARBOXYVINYLTRANSFERASE
13194 (←ここの値を`protein_processed.gro`　に`ligand_fix_GMX.gro`　の値を加えたもの変更)
  312THR      N    1   2.239   3.120  -2.155
  312THR     H1    2   2.167   3.093  -2.219
  312THR     H2    3   2.309   3.173  -2.204
  312THR     H3    4   2.280   3.039  -2.115
  312THR     CA    5   2.181   3.204  -2.049
  312THR     HA    6   2.113   3.147  -2.003
  312THR     CB    7   2.291   3.246  -1.949
  312THR     HB    8   2.350   3.312  -1.995
  312THR    CG2    9   2.242   3.308  -1.818
  312THR   HG21   10   2.321   3.331  -1.761
  312THR   HG22   11   2.189   3.390  -1.837
  312THR   HG23   12   2.185   3.242  -1.769
  312THR    OG1   13   2.364   3.131  -1.914
  312THR    HG1   14   2.436   3.155  -1.850
  312THR      C   15   2.114   3.328  -2.109
  312THR      O   16   2.141   3.361  -2.224
  313GLN      N   17   2.030   3.395  -2.030
  313GLN      H   18   1.996   3.345  -1.950
  313GLN     CA   19   1.984   3.532  -2.047
  313GLN     HA   20   2.051   3.583  -2.101
  313GLN     CB   21   1.850   3.528  -2.122
  313GLN    HB1   22   1.788   3.468  -2.071
  313GLN    HB2   23   1.866   3.489  -2.212
  313GLN     CG   24   1.780   3.663  -2.140
  313GLN    HG1   25   1.838   3.724  -2.195
  313GLN    HG2   26   1.762   3.705  -2.051
  313GLN     CD   27   1.647   3.644  -2.212
  313GLN    OE1   28   1.540   3.644  -2.152
  313GLN    NE2   29   1.651   3.613  -2.339
  313GLN   HE21   30   1.739   3.602  -2.386
  313GLN   HE22   31   1.566   3.599  -2.390
  （中略）
  754VAL      C13135   2.540   3.968  -3.680
  754VAL      O13136   2.545   4.034  -3.576
  755SER      N13137   2.428   3.946  -3.747
  755SER      H13138   2.432   3.882  -3.824
  755SER     CA13139   2.298   4.008  -3.718
  755SER     HA13140   2.308   4.046  -3.626
  755SER     CB13141   2.267   4.118  -3.820
  755SER    HB113142   2.179   4.161  -3.800
  755SER    HB213143   2.266   4.081  -3.913
  755SER     OG13144   2.364   4.221  -3.818
  755SER     HG13145   2.340   4.291  -3.886
  755SER      C13146   2.185   3.905  -3.718
  755SER      O13147   2.203   3.790  -3.760
  756GLN      N13148   2.068   3.948  -3.671
  756GLN      H13149   2.067   4.039  -3.630
  756GLN     CA13150   1.941   3.875  -3.674
  756GLN     HA13151   1.954   3.786  -3.718
  756GLN     CB13152   1.896   3.849  -3.528
  756GLN    HB113153   1.806   3.807  -3.530
  756GLN    HB213154   1.892   3.936  -3.480
  756GLN     CG13155   1.993   3.755  -3.454
  756GLN    HG113156   2.087   3.785  -3.471
  756GLN    HG213157   1.981   3.662  -3.488
  756GLN     CD13158   1.971   3.753  -3.302
  756GLN    OE113159   1.990   3.850  -3.231
  756GLN    NE213160   1.945   3.638  -3.243
  756GLN   HE2113161   1.940   3.554  -3.297
  756GLN   HE2213162   1.930   3.635  -3.144
  756GLN      C13163   1.842   3.956  -3.758
  756GLN    OC113164   1.849   3.934  -3.881
  756GLN    OC213165   1.808   4.066  -3.709
    1  UNL    C    1   0.698   2.601  -2.127(←ここから)
    1  UNL   C1    2   0.758   2.551  -2.013
    1  UNL   C2    3   0.691   2.556  -1.891
    1  UNL   C3    4   0.563   2.611  -1.881
    1  UNL   C4    5   0.500   2.661  -1.999
    1  UNL   C5    6   0.570   2.656  -2.122
    1  UNL   C6    7   0.504   2.615  -1.745
    1  UNL   C7    8   0.388   2.669  -1.704
    1  UNL   C8    9   0.285   2.729  -1.794
    1  UNL    O   10   0.346   2.807  -1.892
    1  UNL   C9   11   0.361   2.718  -1.998
    1  UNL  C10   12   0.183   2.813  -1.716
    1  UNL  C11   13   0.246   2.938  -1.657
    1  UNL   O1   14   0.130   2.737  -1.612
    1  UNL    H   15   0.130   2.643  -1.637
    1  UNL   O2   16   0.513   2.703  -2.234
    1  UNL   H1   17   0.421   2.694  -2.248
    1  UNL   H2   18   0.750   2.598  -2.221
    1  UNL   H3   19   0.856   2.508  -2.018
    1  UNL   H4   20   0.739   2.518  -1.803
    1  UNL   H5   21   0.563   2.569  -1.668
    1  UNL   H6   22   0.368   2.669  -1.599
    1  UNL   H7   23   0.232   2.649  -1.841
    1  UNL   H8   24   0.344   2.769  -2.090
    1  UNL   H9   25   0.292   2.638  -1.986
    1  UNL  H10   26   0.106   2.843  -1.784
    1  UNL  H11   27   0.323   2.911  -1.588
    1  UNL  H12   28   0.288   2.997  -1.736
    1  UNL  H13   29   0.171   2.995  -1.606(←ここまではligand_fix_GMX.groから取ってくる)
   7.38930  10.33800  11.32040
```

topol.topの変更もしてください。

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

; Include chain topologies
#include "topol_Protein_chain_A.itp"
#include "topol_Protein_chain_B.itp"

; Include ligand topology (←ここを追加)
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
UNL		            1 (←ここを追加)
```

これにてタンパク質とリガンドの複合体ファイルの完成です。
