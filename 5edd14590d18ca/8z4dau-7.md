---
title: "📰 NVT（定温）とNPT（定圧定温）の平衡化"
free: false
---

NVT（定温）とNPT（定圧定温）の平衡化を行います。まずは温度を310 Kに徐々に上げる必要があります。

**NVT平衡化（20 ps）**

以下のようにnvt.mdpを設定し、現在のディレクトリに入れてください。**integrator、nsteps、dt、tc-grps、ref_t、gen_temp**を論文用に変更しています。

```
title                   = Protein-ligand complex NVT equilibration
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = sd        (←ここを変更）; Langevin dynamics (stochastic dynamics)
nsteps                  = 10000     (←ここを変更）; 2 * 10000 = 20 ps
dt                      = 0.002     (←ここを変更）; 2 fs
; Output control
nstenergy               = 500   ; save energies every 1.0 ps
nstlog                  = 500   ; update log file every 1.0 ps
nstxout-compressed      = 500   ; save coordinates every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = System    (←ここを変更）; two coupling groups - more accurate
tau_t                   = 0.1                     ; time constant, in ps
ref_t                   = 310   (←ここを変更）; reference temperature, changed to 310 K for each group
; Pressure coupling
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 310   (←ここを変更） ; temperature for Maxwell distribution changed to 310 K
gen_seed                = -1        ; generate a random seed

```

以下を実行してください。

```powershell
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
```

`gmx grompp` コマンドは、シミュレーションの入力ファイル（`.tpr`）を作成するために使用されます。このコマンドでは、エネルギー最小化（`em.gro`）から得られた構造を元に、NVT平衡化の設定ファイル（`nvt.mdp`）を読み込んで、必要なすべての情報を含むバイナリ形式の入力ファイル（`nvt.tpr`）を作成します。

以下を実行して、温度の平衡化を行ってください。

```powershell
gmx mdrun -deffnm nvt
```

**NPT平衡化（100 ps）**

続いて、圧力を平衡化するために、以下のnpt.mdpを作成し、現在のディレクトリに入れてください。**integrator、nsteps、dt、tcoupl、tc-grps、tau_t、ref_t、pcoupl**を論文用に変更しています。

**※pcouplの項目において、チュートリアルではBerendsen  法が設定されていますが、GROMACS2021以降ではC-rescale が推奨されているので、こちらを使用します。**

npt.mdp

```
title                   = Protein-ligand complex NPT equilibration
define                  = -DPOSRES  ; position restrain the protein and ligand
; Run parameters
integrator              = md        (←ここを変更）; leap-frog integrator
nsteps                  = 50000     (←ここを変更）; 2 * 50000 = 100 ps
dt                      = 0.002     (←ここを変更）; 2 fs
; Output control
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
nstxout-compressed      = 500       ; save coordinates every 1.0 ps
; Bond parameters
continuation            = yes       ; continuing from NVT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     (←ここを変更）; modified Berendsen thermostat
tc-grps                 = System   (←ここを変更） ; two coupling groups - more accurate
tau_t                   = 0.1                      (←ここを変更）; time constant, in ps
ref_t                   = 310                      (←ここを変更）; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = C-rescale    (←ここを変更）; 圧力カップリングにC-rescale法を使用
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no
; Velocity generation
gen_vel                 = no        ; velocity generation off after NVT
```

以下を実行してください。

```powershell
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
```

このコマンドは、GROMACSでNPT平衡化（定圧・定温アンサンブル）を実行するための準備を行う`grompp`コマンドです。具体的には、NPTシミュレーションの実行に必要な設定ファイル（`npt.tpr`）を生成します。

続いて、以下を実行し、圧力の平衡化を行います。

```powershell
gmx mdrun -deffnm npt -v
```
