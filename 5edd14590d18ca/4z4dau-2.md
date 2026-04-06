---
title: "📰 タンパク質の準備"
free: false
---
## タンパク質の準備

論文のMaterials and methodsの「Retrieval EPSP synthase crystal structure and preparation for docking studies」の部分フォローしていきます。

ここではタンパク質の準備として、PDBからの抽出やその後のエネルギー最小化を行います。
PDBは化合物との複合体になっていることが多いので、少なからずタンパク質の構造が変化しています。

そのため、エネルギー最小化でタンパク質を化合物が結合していない状態に戻しています。
（結合状態を再現した状態で、in silico創薬していきたい場合などはエネルギー最小化を行わない場合もあります。）

まずは[5BUFのPDB](https://www.rcsb.org/structure/5BUF)からPDBファイルをダウンロードしてください。

![](https://storage.googleapis.com/zenn-user-upload/9e3c91ed4990-20241016.jpg)


※タンパク質の構造情報がない方は
**第2章：Homology modelingによる構造予測**
から予測した構造を使うと良いと思います。

もしくはAlphaFoldやBoltzなどで予測した構造を使用し、試してみてください。

ここではOpenMMを使った方法で行います。

## OpenMMとは

OpenMM（オープンエムエム）は、Python や C++ から呼び出せる **GPU対応の分子動力学（MD）シミュレーション・ライブラリ**です。CUDA または OpenCL を活用でき、数万原子規模の系でも実用的な速度でカスタム計算を行えます。

主な特長:

- **Python API** を介して、`ForceField` や `Simulation`、カスタム積分器などを柔軟に組み立て可能です。
- **GPU 有効化** により、CUDA／OpenCL ベースの高速化が進んでおり、複数GPUでの並列実行もサポートされます。
- **構造最適化（minimize）** やサンプリング、自由エネルギー計算などの基本機能を備えています。

ライセンスは、**API や CPU リファレンス部分が MIT ライセンス**（商用・改変可）で公開されており、CUDA／OpenCL プラットフォーム部分は **LGPL ライセンス** です。

**総じて、OpenMM は、柔軟なカスタマイズ性とGPUによる高性能を両立し、研究用途から創薬シミュレーションまで幅広く利用できる強力なMDツールキット**です。

ここではOpenMMを使って、エネルギーを最小化してみましょう！

## 環境構築

可視化ツールであるPyMOL（オープンソース版）と共に、OpenMMの環境構築をしましょう！
以下をターミナルで実行してください。

```python
conda create -n protein_preparation python=3.10 -y
conda activate protein_preparation

conda install -c conda-forge pymol-open-source
conda install -c conda-forge openmm

pymol
```

### 詳細コード解説

**仮想環境の作成と起動**

```bash
conda create -n protein_preparation python=3.10 -y
```

- **仮想環境 `protein_preparation` をPython 3.10で作成**
- `y`：確認なしで自動的にインストールを進める

```bash
conda activate protein_preparation
```

- 作成した仮想環境を有効化
- これ以降の操作はこの環境内で行われます

---

**PyMOL オープンソース版のインストール**

```bash
conda install -c conda-forge pymol-open-source
```

- conda-forge チャンネルから**PyMOLオープンソース版**をインストール
- 商用ライセンス不要で、拡張機能やスクリプト利用も可

**注意点：**

- Schrödinger版とは異なり、GUIの見た目や機能が少しシンプルです
- Macでは起動時に`XQuartz`（X11）が必要になることがあります

---

**OpenMMのインストール**

```bash
conda install -c conda-forge openmm
```

- 分子力学シミュレーションエンジン OpenMM をインストール
- Pythonから呼び出して**エネルギー最小化やMD**を行うことができます

---

**PyMOLの起動**

```bash
pymol
```

- PyMOLを起動します（GUIが立ち上がります）
- PyMOLの下部ターミナルでPythonコードやスクリプトを実行できます

## PyMOLを使ったタンパク質の準備

---

まずPyMOLが開かれると、以下の画面になっていると思います。左上のFIleからタンパク質をロードしてください。

![](https://storage.googleapis.com/zenn-user-upload/860ce7b4262b-20250809.png)

まずは水分子を消しましょう。

![](https://storage.googleapis.com/zenn-user-upload/d57081e5887f-20250809.png)

次に余分なイオンとリガンドを消していきます。

上タブのDisplay→Sequenceを押して、配列を出しましょう。

![](https://storage.googleapis.com/zenn-user-upload/41bf7dd0e29b-20250809.png)

イオンとリガンドを選択し、右側のremove atomsから消してください。

![](https://storage.googleapis.com/zenn-user-upload/c201d256337e-20250809.png)

左上のFile→Export Moleculesを押し、5buf_no_water_ligand.pdb　として保存してください。

## OpenMMを使ったエネルギー最小化

---

続いて、OpenMMを使って、エネルギーを最小化していきます。
第9章のAnacondaによる環境構築を参考にしながら、Python環境を整えてください。
その後以下のコードを実行してみてください。

全コードはこちら

※ `/path/to/` の部分はご自身のファイルパスに置き換えてください。

```python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# 1. PDB読み込み
pdb = PDBFile('/path/to/5buf_no_water_ligand.pdb')

# 2. ForceField定義
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

# 3. Modellerを使って水素原子を追加
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

# 4. システム構築
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)

# 5. 積分器の定義
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# 6. シミュレーション準備
simulation = Simulation(modeller.topology, system, integrator)

# 7. 位置を設定
simulation.context.setPositions(modeller.positions)

# 8. エネルギー最小化を実行
print("エネルギー最小化を実行中...")
simulation.minimizeEnergy()

# もっとしたい場合はこちら
#simulation.minimizeEnergy(
#    tolerance=10 * kilojoule_per_mole / nanometer,
#    maxIterations=1000
#)
print("エネルギー最小化完了！")

# 9. 最小化後の構造を取得
positions = simulation.context.getState(getPositions=True).getPositions()

# 10. PDBファイルに出力
with open('minimized.pdb', 'w') as output_file:
    PDBFile.writeFile(modeller.topology, positions, output_file)

print("最小化された構造を 'minimized.pdb' に保存しました。")
```

### 詳細コード

---

**必要なモジュールのインポート**

```python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
```

- `openmm.app`：構造データの読み込み・操作、シミュレーション設定に使います。
- `openmm`：シミュレーションエンジンそのものです。
- `openmm.unit`：温度や長さなどに単位（K, nmなど）を付けるために使います。
- `sys.stdout`：ログやメッセージの表示先に使うこともあります（今回は未使用でもOK）。

---

**PDBファイルの読み込み**

```python
pdb = PDBFile('/path/to/5buf_no_water_ligand.pdb')
```

- `PDBFile`で構造ファイル（`.pdb`）を読み込みます。
- ここでは水やリガンドが除かれたPDBを使用。

---

**ForceFieldの定義**

```python
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
```

- `amber99sb.xml`：タンパク質に使われる力場。
- `tip3p.xml`：水分子に対する力場（今回はまだ水は加えていませんが、後で使う可能性あり）。

---

**Modellerを使って水素原子を追加**

```python
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
```

- PDB構造には通常**水素が含まれていない**ため、自動で水素を追加します。
- `Modeller`オブジェクトで、構造の修正を行えます。

---

**システムの構築**

```python
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)
```

- 力場とトポロジー情報から力の計算モデル（System）を生成。
- `PME`：長距離の電荷相互作用を正確に扱う方法。
- `HBonds`：水素結合のみを固定（計算の高速化と安定化）。

---

**積分器（Integrator）の設定**

```python
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
```

- 温度を制御しながら時間発展を行う**Langevin Dynamics**を使用。
- 300Kで、1/ps の摩擦係数、2 fs のタイムステップ。

---

**Simulationオブジェクトの作成**

```python
simulation = Simulation(modeller.topology, system, integrator)
```

- 上記で作成したトポロジー・力場・積分器を使って、シミュレーション本体を準備。

---

**初期構造の設定**

```python
simulation.context.setPositions(modeller.positions)
```

- 初期の原子座標をシミュレーションに反映させます。

---

**エネルギー最小化を実行**

```python
print("エネルギー最小化を実行中...")
simulation.minimizeEnergy()

#もっと厳密にエネルギー最小化をしたい場合はこちら
#simulation.minimizeEnergy(
#    tolerance=10 * kilojoule_per_mole / nanometer,
#    maxIterations=1000
#)
print("エネルギー最小化完了！")

```

- 原子間の衝突や不自然な構造を修正（エネルギー最小化）します。
- 最小化によって、より現実的な構造へと近づけます。

パラメータを調整することで、さらに厳密にエネルギー最小化を行えます。

 `maxIterations`　の値を変更してみてください。

---

**最小化後の構造を取得**

```python
positions = simulation.context.getState(getPositions=True).getPositions()
```

- 最小化された座標を取得します。
- この座標を使って新しいPDBファイルを書き出せます。

---

**最小化構造のPDBファイル出力**

```python
with open('minimized.pdb', 'w') as output_file:
    PDBFile.writeFile(modeller.topology, positions, output_file)

print("最小化された構造を 'minimized.pdb' に保存しました。")
```

- `writeFile()` でPDB形式として保存します。
- これにより、**水素付き＆最小化済みの構造が `minimized.pdb` に出力されます**。
