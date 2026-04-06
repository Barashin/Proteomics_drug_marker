---
title: "📰 エネルギー最小化"
free: false
---

続いて、プロトコルでは1000ステップのエネルギー最小化を行います。

以下をテキストファイルでコピーして、em.mdpファイルを作成してください。

nstepsの部分を1000ステップに変更します。

```
; LINES STARTING WITH ';' ARE COMMENTS
title		    = Minimization	; Title of run

; Parameters describing what to do, when to stop and what to save
integrator	    = steep		; Algorithm (steep = steepest descent minimization)
emtol		    = 1000.0  	; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps		    = 1000	(←ここを変更)  	; Maximum number of (minimization) steps to perform 

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		    ; Method to determine neighbor list (simple, grid)
rlist		    = 1.2		    ; Cut-off for making neighbor list (short range forces)
coulombtype	    = PME		    ; Treatment of long range electrostatic interactions
rcoulomb	    = 1.2		    ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw		    = 1.2		    ; long range Van der Waals cut-off
pbc             = xyz 		    ; Periodic Boundary Conditions
DispCorr        = no
```

`em.mdp`ファイルを現在のディレクトリに置いて、エネルギー最小化を行います。

```powershell
gmx grompp -f em.mdp -c solvated_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

`gmx grompp -f em.mdp -c solvated_ions.gro -p topol.top -o em.tpr`

- `gmx grompp`: GROMACSの前処理コマンドです。力場のパラメータファイルや構造ファイルを使って、シミュレーションを実行するための入力ファイル（`.tpr`ファイル）を作成します。
- `f em.mdp`: シミュレーションのパラメータファイル（この場合、エネルギー最小化の設定が記載された`em.mdp`ファイル）を指定します。
- `c solvated_ions.gro`: エネルギー最小化の初期構造として使用するGROMACS形式の座標ファイル（溶媒とイオンが含まれた`solvated_ions.gro`）。
- `p topol.top`: トポロジーファイル（シミュレーションに使用する分子の構造や結合情報などが記載されたファイル）を指定します。
- `o em.tpr`: gromppが作成する出力ファイル。`em.tpr`は、GROMACSがシミュレーションを実行するためのバイナリ形式のファイルです。

`gmx mdrun -v -deffnm em`

- `gmx mdrun`: 実際にシミュレーションを実行するためのコマンドです。`grompp`で作成された`.tpr`ファイルを使用して、エネルギー最小化やMDシミュレーションを実行します。
- `v`: 実行中に詳細な出力を表示するオプションです。進行状況をリアルタイムで確認できます。
- `deffnm em`: 出力ファイルの接頭辞（デフォルトファイル名）を指定します。この場合、すべての出力ファイル（エネルギー最小化のログ、座標、エネルギーファイルなど）が`em`という名前で生成されます（例：`em.log`, `em.gro`, `em.edr`など）

以下のように出力できればOKです。

```powershell

Steepest Descents did not converge to Fmax < 1000 in 1001 steps.
Potential Energy  = -2.6026952e+06
Maximum force     =  8.7692803e+03 on atom 6580
Norm of force     =  4.0724692e+01

GROMACS reminds you: "Everybody's Good Enough For Some Change" (LIVE)
```

Steepest Descents did not converge to Fmax ＜ １０００ 

**というメッセージがあるため、GROMACSのエネルギー最小化過程で設定した停止条件に達しなかったことを示しています。**

そのため、ステップ数をさらに上げるべきかもしれませんが、今回は論文にそって、1000ステップのままでいきます。

