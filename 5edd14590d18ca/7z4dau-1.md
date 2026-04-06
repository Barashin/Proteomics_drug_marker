---
title: "📰 WindowsへのLinuxの導入"
free: false
---

コントロールパネル→プログラム→Windowsの機能の有効化または無効化を開いてください。

![](https://storage.googleapis.com/zenn-user-upload/466474044d81-20241116.png)

以下の`Linux用Windowsサブシステム`にチェックが入っていることを確認してください。

![](https://storage.googleapis.com/zenn-user-upload/d9edaa219792-20241116.jpg)

続いて、Ubuntuのアプリをダウンロードしてください。

![](https://storage.googleapis.com/zenn-user-upload/3c6fda1bee03-20241116.jpg)

Ubuntuを開いて、以下のようにusernameとpasswordを設定してください。

![](https://storage.googleapis.com/zenn-user-upload/8402bad65ae8-20241116.png)

これが成功すると、ubuntuが使えるようになります！

以下をubuntuのCUIで実行してください。

そうすると、windowsシステムのデスクトップが作業ディレクトリになります。

```powershell
cd /mnt/c/users/(あなたのユーザー名）/desktop
```

![](https://storage.googleapis.com/zenn-user-upload/f6b52268a3d3-20241116.png)

`/mnt`は、LinuxやUnix系のオペレーティングシステムで、外部のファイルシステム（例えば、外付けハードドライブやUSBドライブなど）を一時的にマウントするための標準的なディレクトリです。`/mnt`を使うことで、これらの外部デバイスにアクセスし、ファイルの読み書きができるようになります。
