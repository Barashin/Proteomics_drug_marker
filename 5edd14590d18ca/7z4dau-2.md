---
title: "📰 GROMACSのインストール"
free: false
---

GROMACSのインストール方法を以下で説明します。

Micromambaを使って仮想環境を構築し、GROMACSをインストールします。

Micromambaは、Pythonや各種ソフトウェアの仮想環境を管理するツールです。
condaの機能を持ちながら、より高速で軽量なのが特徴です。


```bash
# Micromambaをインストールするディレクトリを作成
mkdir -p ~/micromamba_env

# Micromambaをダウンロードしてインストール（-b: バッチモード, -p: インストール先指定）
curl -Ls https://micro.mamba.pm/install.sh | bash -s -- -b -p ~/micromamba_env

# .bashrcファイルを再読み込みしてPATHを更新
source ~/.bashrc

# GROMACS用の仮想環境を作成(CPUのみ、もしくはMacでGPUありOpenClの方)
micromamba create -n gromacs_env -c conda-forge "gromacs=2024" openbabel -y

# GROMACS用の仮想環境を作成(WindowsでGPU:CUDAを使う場合)
micromamba create -n gromacs_env -c conda-forge "gromacs=2024.*=*cuda*" openbabel -y
echo 'alias gmx=gmx_mpi' >> ~/.bashrc
source ~/.bashrc

# 作成した環境をアクティベート
micromamba activate gromacs_env

#gromaxのインストール確認(GPUを確認してください）
gmx pdb2gmx --version

#acpype,ambertools のインストール
pip install acpype biopython
micromamba install -c conda-forge ambertools -y

```

### コード詳細

`mkdir -p ~/micromamba_env`

これは、`micromamba_env`という名前のディレクトリ（フォルダ）を、ユーザーのホームディレクトリ（`~/`）に作成するコマンドです。`-p`オプションは、もし親ディレクトリが存在しない場合でも、それらも一緒に作成することを意味します。このディレクトリは、後で`micromamba`の実行ファイルをインストールするために使われます。

---

`curl -Ls https://micro.mamba.pm/install.sh | bash -s -- -b -p ~/micromamba_env`

このコマンドは、`micromamba`をインストールするために実行されます。

- `curl -Ls https://micro.mamba.pm/install.sh`：指定されたURLから`micromamba`のインストールスクリプトをダウンロードします。`L`はリダイレクトを追跡し、`s`は進行状況表示を抑制します。
- `|`：パイプ演算子で、前のコマンドの出力を次のコマンドの入力に渡します。
- `bash -s -- -b -p ~/micromamba_env`：ダウンロードしたスクリプトを`bash`で実行します。
    - `s`：スクリプトが標準入力から読み込まれることを示します。
    - `-`：後続の引数がスクリプトに渡されることを示します。
    - `b`：バッチモードでインストールを実行し、ユーザーとの対話を省略します。
    - `p ~/micromamba_env`：インストール先を`~/micromamba_env`ディレクトリに指定します。

---

`source ~/.bashrc`

これは、`.bashrc`ファイルに書かれている設定を現在のシェルセッションに読み込むコマンドです。`micromamba`のインストール時に、実行ファイルへのパス（`PATH`）が`.bashrc`に追加されます。このコマンドを実行することで、新しく追加された`micromamba`コマンドがすぐに使えるようになります。

---

`micromamba create -n gromacs_env -c conda-forge gromacs openbabel -y`

このコマンドは、GROMACSとOpen Babelがインストールされた、`gromacs_env`という名前の新しい仮想環境を作成します。

- `micromamba create`：新しい仮想環境を作成するコマンドです。
- `n gromacs_env`：作成する仮想環境の名前を`gromacs_env`と指定します。
- `c conda-forge`：パッケージを`conda-forge`というチャネル（リポジトリ）から取得するように指定します。`conda-forge`は、科学計算用のパッケージが豊富に揃ったコミュニティ主導のチャネルです。
- `"gromacs=2024" openbabel`：インストールするパッケージの名前です。
- `y`：確認メッセージを省略し、自動的に処理を進めます。

---

`micromamba create -n gromacs_env -c conda-forge "gromacs=2024.*=*cuda*" openbabel -y`

こちらは、GPU（CUDA）を利用するためのGROMACSをインストールするためのコマンドです。

- `"gromacs=2024.*=*cuda*"`：GROMACSのバージョンが2024系で、かつCUDAをサポートするビルドをインストールすることを指定します。

`echo 'alias gmx=gmx_mpi' >> ~/.bashrc`

このコマンドは、`gmx`コマンドを実行した際に、自動的に`gmx_mpi`が呼び出されるように、エイリアスを設定します。

- `echo 'alias gmx=gmx_mpi'`：`alias gmx=gmx_mpi`という文字列を標準出力に出力します。
- `>> ~/.bashrc`：その出力を`.bashrc`ファイルの末尾に追記します。
これにより、並列計算に特化した`gmx_mpi`がデフォルトで使われるようになります。

`source ~/.bashrc`

エイリアスをすぐに有効にするために、`.bashrc`ファイルを再読み込みします。

---

`micromamba activate gromacs_env`

これは、先ほど作成した`gromacs_env`という仮想環境をアクティベート（有効化）するコマンドです。このコマンドを実行すると、ターミナルのプロンプトが変わり、この仮想環境にインストールされたパッケージやコマンドが使える状態になります。

---

`gmx pdb2gmx --version`

GROMACSが正しくインストールされているか、またGPUが認識されているかを確認するためのコマンドです。`gmx pdb2gmx`はGROMACSのサブコマンドの一つで、`--version`オプションを付けると、GROMACSのバージョン情報と、OpenMPやCUDAなどの並列化オプションが有効になっているかどうかが表示されます。

---

`pip install acpype biopython`

`acpype`と`biopython`という2つのPythonパッケージをインストールします。

- `acpype`：分子の力場ファイルを生成するのに役立つツールです。
- `biopython`：バイオインフォマティクス分野で広く使われるPythonライブラリです。
これらのパッケージは、仮想環境がアクティベートされている状態で実行すると、その仮想環境内にインストールされます。

---

`micromamba install -c conda-forge ambertools -y`

これは、`micromamba`を使って`ambertools`パッケージをインストールするコマンドです。

- `ambertools`：分子動力学シミュレーションなどで使われる、AMBER力場のツールキットです。
`pip`ではなく`micromamba`（または`conda`）でインストールすることで、依存関係を適切に解決し、環境の整合性を保つことができます。
