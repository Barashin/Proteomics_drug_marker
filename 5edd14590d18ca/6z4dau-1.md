---
title: "📰 PyRxによる化合物ライブラリのエネルギー最小化"
free: false
---

※本章で使用するPyRxですが、Macでは使用できないようです。Macの場合は
**第7章：Sminaを使ったin silicoスクリーニング**
をご参照ください。

## PyRxのinstall

まず[PyRXのダウンロードページ](https://pyrx.sourceforge.io/downloads)からwindows版のFree 0.8をダウンロードします。

![](https://storage.googleapis.com/zenn-user-upload/1785d4c0a04a-20241105.png)

ダウンロードしたPyRx-0.8-Setup.exeからPyRxをインストールしてください。特に特別なことはせず、Nextを押してインストールしてください。

## PyRx内での低分子のエネルギー最小化

まずは忘れないうちにoutputディレクトリを決めておきます。

PyRxを開き、Edit→Preferencesを押してください。

![](https://storage.googleapis.com/zenn-user-upload/a581947c5749-20241105.jpg)

以下のworkspaceに任意のディレクトリを設定してください。これがOutputディレクトリとなります。

![](https://storage.googleapis.com/zenn-user-upload/a8040a264861-20241105.jpg)

続いて、フィルタリングしたライブラリを取り込んでいきます。

File→Importを押してください。

![](https://storage.googleapis.com/zenn-user-upload/898d3f14f95d-20241105.jpg)

SDFを押して、フィルタリングした`CMNPD-filtered` をアップロードしてください。
（※第5章で構築したライブラリとは異なりますが、同じ手順でできますので、そのまま続けてください）

![](https://storage.googleapis.com/zenn-user-upload/cf4fe19eb942-20241105.png)

以下のようにアップロードできたと思うので、すべての化合物が入っているか確認してください。

`Open Babel`のタブのところにあります。

![](https://storage.googleapis.com/zenn-user-upload/a3e9b8f99673-20241105.jpg)

適当な化合物に対して、右クリックを押し、すべての低分子化合物のエネルギーを最小化します。

論文にはエネルギー最小化の方法が書かれていないので、デフォルトのuffで行います。

![](https://storage.googleapis.com/zenn-user-upload/f95d9b779607-20241105.jpg)

Force Fieldの変更の仕方

他のForce Fieldで行いたい方は、以下の画像の右上にあるMを押して、一つを別の方法でエネルギー最小化を行った後、再度Minimize Allを行ってください。

![](https://storage.googleapis.com/zenn-user-upload/d0b6c0f07842-20241105.jpg)

Estimated timeは15分くらいですが、だいたい30分くらいかかります。

![](https://storage.googleapis.com/zenn-user-upload/5852ddaaa339-20241105.png)

続いて、すべてのファイルをpdbqtファイルに変換します。

Estimated timeは15分くらいですが、だいたい１時間くらいかかります。

![](https://storage.googleapis.com/zenn-user-upload/e920c1a9ee2e-20241105.png)

![](https://storage.googleapis.com/zenn-user-upload/71c9b0dbfe01-20241105.png)

これで化合物ライブラリの準備は終了しました！
