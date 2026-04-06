---
title: "📰 MacでのGROMACSのインストール"
free: false
---
(こちらは前Chapterでうまくいかなかった時に試してみてください。)

## CPUのみを用いる場合の簡単なGROMACSのインストール


GROMACSを使うにあたり、[Homebrew](https://brew.sh/index_ja)をまずターミナルからインストールします。

Homebrewは、macOS向けのパッケージ管理システムです。パッケージ管理システムは、ソフトウェアのインストールやアップデートを簡単に行えるツールであり、開発者やユーザーにとって便利です。

Homebrewを使用すると、ターミナルを介してコマンドを実行することで、様々なソフトウェアパッケージをインストールできます。これにより、開発ツール、プログラミング言語、データベース、ライブラリなど、さまざまなソフトウェアを簡単に管理できます。

Homebrewを使用すると、Mac上でGROMACSをインストールすることができます。Homebrewを使ってGROMACSをインストールすると、依存関係の解決やバージョン管理などが自動的に行われます。

まずターミナルで

```sql
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

と打ち、Homebrewをインストールします。

続いて

```sql
brew install gromacs
```

と打ち、gromacsをインストールします。

```sql
gmx
```

と打ち、以下のような記述が返ってきます。

```sql
:-) GROMACS - gmx, 2023.1-Homebrew (-:

Executable:   /opt/homebrew/bin/../Cellar/gromacs/2023.1/bin/gmx
Data prefix:  /opt/homebrew/bin/../Cellar/gromacs/2023.1
Working dir:  /Users/kshinba
Command line:
  gmx

SYNOPSIS

gmx [-[no]h] [-[no]quiet] [-[no]version] [-[no]copyright] [-nice <int>]
    [-[no]backup]

OPTIONS

Other options:

 -[no]h                     (no)
           Print help and quit
 -[no]quiet                 (no)
           Do not print common startup info or quotes
 -[no]version               (no)
           Print extended version information and quit
 -[no]copyright             (no)
           Print copyright information on startup
 -nice   <int>              (19)
           Set the nicelevel (default depends on command)
 -[no]backup                (yes)
           Write backups if output files exist

Additional help is available on the following topics:
    commands    List of available commands
    selections  Selection syntax and usage
To access the help, use 'gmx help <topic>'.
For help on a command, use 'gmx help <command>'.

GROMACS reminds you: "Correctomundo" (Pulp Fiction)
```

これが出てくるときちんとHomebrewとGROMACSのインストールが完了しています。

## GPUを用いる場合のGROMACSのインストール

こちらではGPUを使う場合を説明します。GPUがある場合はこちらを行った方が断然計算スピードが速くなります。

### 1. 必要なライブラリのインストール

```
# 必要なパッケージのインストール
sudo apt-get update
sudo apt-get install -y build-essential cmake wget
```

- `cmake`：ビルドツール。GROMACSのコンパイルに使用。
- `fftw`：高速フーリエ変換ライブラリ。分子動力学の計算効率を向上。
- `openmpi`：並列計算ライブラリ。マルチスレッド処理を高速化。

### 2. GROMACSのインストール（GPU版）

```bash
# GROMACS 2025.1のダウンロード
wget ftp://ftp.gromacs.org/gromacs/gromacs-2025.1.tar.gz

# アーカイブの展開
tar -xvzf gromacs-2025.1.tar.gz
cd gromacs-2025.1
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_MPI=OFF -DGMX_OPENMP=ON -DGMX_GPU=OpenCL -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs
make -j$(sysctl -n hw.logicalcpu)
sudo make install
```

- `wget`：GROMACSのソースコードをダウンロード。
- `tar -xvzf`：ダウンロードしたファイルを解凍。
- `mkdir build`：ビルド用のディレクトリを作成。
- `cmake ..`：CMakeを使用してコンパイルオプションを設定。
- `make -j$(sysctl -n hw.logicalcpu)`：利用可能なCPUコア数を自動的に検出し並列コンパイル。
- `sudo make install`：GROMACSをインストール。

GROMACSのビルド時に使用するCMakeオプションについて、それぞれの引数の意味と役割を以下にまとめます。

1. **`DGMX_BUILD_OWN_FFTW=ON`**
    - **説明**: GROMACSは高速フーリエ変換（FFT）を多用します。このオプションを有効にすると、GROMACS自身がFFTWライブラリをソースからダウンロードし、ビルドします。
    - **使用目的**: システムに適切なFFTWライブラリがインストールされていない場合や、GROMACSと最適な互換性を持つFFTWを使用したい場合に有効です。
    - **参考**: [GROMACS 2024.4 インストールガイド](https://manual.gromacs.org/2024.4/install-guide/index.html)
2. **`DGMX_MPI=OFF`**
    - **説明**: GROMACSでMPI（Message Passing Interface）を使用した並列計算を有効にします。
    - **使用目的**: 複数の計算ノード間での並列計算を行いたい場合に設定します。ただし、シングルノード内での並列計算には不要です。なので、今回はOFFにしています。[ArchWiki](https://wiki.archlinux.jp/index.php/GROMACS?utm_source=chatgpt.com)
    - **注意**: このオプションを有効にするには、OpenMPIやMPICHなどのMPIライブラリが事前にインストールされている必要があります。
    - **参考**: [ArchWikiのGROMACSページ](https://wiki.archlinux.jp/index.php/GROMACS)[ArchWiki](https://wiki.archlinux.jp/index.php/GROMACS?utm_source=chatgpt.com)
3. **`DGMX_OPENMP=ON`**
    - **説明**: GROMACSでOpenMPを使用したマルチスレッド並列計算を有効にします。
    - **使用目的**: 単一の計算ノード内で、複数のCPUコアを利用して計算速度を向上させる場合に設定します。
    - **参考**: [GROMACS 2024.4 インストールガイド](https://manual.gromacs.org/2024.4/install-guide/index.html)
4. **`DGMX_GPU=OpenCL`**
    - **説明**: GROMACSでのGPU（Graphics Processing Unit）サポートを有効にします。
    - **使用目的**: GPUを使用環境でのビルド時に設定します。
    - **参考**: [GROMACS 2024.4 インストールガイド](https://manual.gromacs.org/2024.4/install-guide/index.html)
5. **`DCMAKE_INSTALL_PREFIX=/usr/local/gromacs`**
    - **説明**: GROMACSをインストールするディレクトリを指定します。
    - **使用目的**: デフォルト以外の場所にGROMACSをインストールしたい場合に設定します。
    - **補足**: インストール先のディレクトリに書き込み権限が必要です。
    - **参考**: [GROMACSをVisual Studioでビルドする話](https://sage-t.tumblr.com/post/83732368436/gromacs-visual-studio)[Tumblr](https://sage-t.tumblr.com/post/83732368436/gromacs-visual-studio?utm_source=chatgpt.com)
- `make -j$(sysctl -n hw.logicalcpu)`：利用可能なCPUコア数を自動的に検出し並列コンパイル。
- `sudo make install`：GROMACSをインストール。

環境変数を設定します。

```
echo 'source /usr/local/gromacs/bin/GMXRC' >> ~/.bash_profile
source ~/.bash_profile
```

- `source /usr/local/gromacs/bin/GMXRC`：GROMACSの環境変数を適用。

動作確認を行います。以下のコマンドで`GPU support: OpenCL`　と表示されれば、GPUを認識しています。

```
gmx --version
```

## 参考記事

https://qiita.com/Masa0127/items/b0d5e0b27c4a96ae2ae0

https://manual.gromacs.org/2024.4/install-guide/index.html
