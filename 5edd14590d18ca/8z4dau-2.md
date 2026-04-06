---
title: "📰 タンパク質の準備"
free: false
---

続いて、protein.pdbがあるディレクトリに移動して、以下を実行していきます。

GROMACSでタンパク質をパラメータ化するために`pdb2gmx`を使用します。

**タンパク質のパラメータ化**：
GROMACSのFf14SBに相当するAMBER力場を使うため、以下のコマンドを実行します。GROMACSにはAMBER力場が組み込まれています。このコマンドでは、Ff14SBに相当する`amber99sb-ildn`力場を指定し、TIP3P水モデルを使用します。

```powershell
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -ignh
```

`gmx pdb2gmx`

`gmx`はGROMACSのコマンドラインツールのプレフィックスで、`pdb2gmx`はPDBファイルを読み込み、GROMACSで使用できる形式に変換するコマンドです。具体的には、PDBファイルから力場を選択して、適切なパラメータ（結合情報、水素原子の追加、末端処理など）を設定して出力します。

`f protein.pdb`

`f`は「ファイル」を意味するオプションです。後ろに続く`protein.pdb`は入力ファイル名で、ここではPDB形式のタンパク質構造ファイルを指定しています。このファイルには、タンパク質の三次元構造情報が記述されています。 `o protein_processed.gro`

`o`は「出力ファイル」を指定するオプションです。`protein_processed.gro`は出力されるGROMACS用のファイル名で、`.gro`形式の座標ファイルとして保存されます。`.gro`ファイルは、GROMACSのシミュレーションで使われる標準的なフォーマットで、分子の座標とボックスサイズの情報が含まれています。

 `ignh`

このオプションは「入力ファイル内の水素原子を無視する」ことを意味します。PDBファイルにはしばしば不完全な水素原子の情報が含まれていますが、このオプションを使用することで、既存の水素原子を無視し、GROMACSが必要に応じて自動的に水素原子を追加してくれます。

以下のように出力します。

```powershell
                      :-) GROMACS - gmx pdb2gmx, 2024 (-:

Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /home/shizuku/protein_ligand_complex_ver2
Command line:
  gmx pdb2gmx -f protein.pdb -o protein_processed.gro

Select the Force Field:

From '/usr/local/gromacs/share/gromacs/top':

 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)

 2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)

 3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)

 4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)

 5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)

 6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)

 7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)

 8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)

 9: GROMOS96 43a1 force field

10: GROMOS96 43a2 force field (improved alkane dihedrals)

11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)

12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)

13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)

14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)

15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)
```

`amber99sb-ildn`　を指定するので、6を選択します。

```powershell
Using the Amber99sb-ildn force field in directory amber99sb-ildn.ff
Opening force field file /usr/local/gromacs/share/gromacs/top/amber99sb-ildn.ff/watermodels.dat

Select the Water Model:

 1: TIP3P     TIP 3-point, recommended

 2: TIP4P     TIP 4-point

 3: TIP4P-Ew  TIP 4-point optimized with Ewald

 4: TIP5P     TIP 5-point (see https://gitlab.com/gromacs/gromacs/-/issues/1348 for issues)

 5: SPC       simple point charge

 6: SPC/E     extended simple point charge

 7: None
```

`tip3p`　を指定するので、1を選択します。

以下のファイルが出力されればOKです！

![](https://storage.googleapis.com/zenn-user-upload/eb2b1d6f2721-20241125.png)
