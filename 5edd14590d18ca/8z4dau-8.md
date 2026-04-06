---
title: "📰 MDシミュレーション"
free: false
---

生産MDシミュレーションを行います。GPUが設定されていない場合は10日かかるので、状況に応じて変更してください。

md.mdpを以下のように設定してください。**nsteps、tc-grps、ref_t、pcoupl**を[GROMACSのチュートリアル](http://www.mdtutorials.com/gmx/complex/08_MD.html)から変更しています。

md.mdp

```
title                   = Protein-ligand complex MD simulation
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000000  (←ここを変更）; 2 * 50000000 = 100000 ps (100 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save coordinates every 10.0 ps
; Bond parameters
continuation            = yes       ; continuing from NPT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2       ; neighbor search cutoff in nm
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0       ; start of the switching function for vdW
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2       ; cutoff for Coulomb interactions in nm
pme_order               = 4         ; cubic interpolation for PME
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = System    (←ここを変更）; two coupling groups - more accurate
tau_t                   = 0.1                       ; time constant, in ps
ref_t                   = 310           (←ここを変更）; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = C-rescale    (←ここを変更）; 圧力カップリングにC-rescale法を使用
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no
; Velocity generation
gen_vel                 = no        ; continuing from NPT equilibration

```

以下を実行してください

```powershell
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
```

このコマンドは、GROMACSのコマンドラインツールである`grompp`を使用して、MDシミュレーションの入力ファイル（`.tpr`ファイル）を生成するためのものです。

以下でMDシミュレーションを開始してください。

```powershell
gmx mdrun -deffnm md_0_10 -v
```

終了すると、以下のような出力が出ます。

```powershell
step 50000000, remaining wall clock time:     0 s
               Core t (s)   Wall t (s)        (%)
       Time:  1762898.003    73454.087     2400.0
                         20h24:14
                 (ns/day)    (hour/ns)
Performance:      117.624        0.204

GROMACS reminds you: "I don't fear death because I don't fear anything I don't understand." (Hedy Lamarr)
```

もし途中でシミュレーションが終了していた場合は、以下のコードで再開できます。

```powershell
gmx mdrun -s md_0_10.tpr -cpi md_0_10.cpt -deffnm md_0_10
```
