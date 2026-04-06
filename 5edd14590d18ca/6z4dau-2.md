---
title: "📰 In silicoスクリーニング"
free: false
---

用意したタンパク質とエネルギー最小化した化合物ライブラリを設定します。

Vina WizardタブでStartを教えてください

![](https://storage.googleapis.com/zenn-user-upload/28bb215ebfce-20241105.jpg)

Add Macromoleculeからエネルギー最小化したタンパク質をアップロードします。

![](https://storage.googleapis.com/zenn-user-upload/21929f1fdc0e-20241105.png)

アップロードするとタンパク質が表示されます。

![](https://storage.googleapis.com/zenn-user-upload/1237206e5071-20241105.png)

アップロードしたタンパク質をクリックします。

![](https://storage.googleapis.com/zenn-user-upload/481c615a7f05-20241105.png)

そうすると、`Molecule`タブで各アミノ酸が設定されるので、Ser550 OG (論文ではSer238 OG、アミノ酸開始が１なので。)をクリックして、マークしてください。後で、ドッキング領域を決めるGrid Boxの目安となります。

![](https://storage.googleapis.com/zenn-user-upload/be6cab43d1c2-20241105.jpg)

その後、すべてリガンドを選択して、`Forward`を押してください。

リガンドを選択したとき、リガンドの数が少し減っているかと思います。

これはpdbqtに変換したときにエラーしたものが、含まれていないためです。

スクリーニングの段階なので、とりあえずは無視して先に進みます。

![](https://storage.googleapis.com/zenn-user-upload/c7305316c09c-20241105.jpg)

続いて、Grid Boxを決めていきます。

だいたい論文で定義してある値に揃えてください。

先ほど印をつけたSer550のOGアトムが真ん中に来るようにGrid Boxを設定してください。

X-axis: −0.544 Å, Y-axis: 28.321 Å, and Z-axis: −12.451 Åになるべく合うようにCenterを微調整してください。

以下の図の左上に座標のアイコンがあるので、微調整はそれを参考にすると少しわかりやすいです。Angstromのｘ、ｙ、ｚは１５くらいに設定してください。

![](https://storage.googleapis.com/zenn-user-upload/a12a4e82dba2-20241105.jpg)

論文には100 iterationという記述があり、これが少し不明です。[他の論文を読む](https://app.paperpile.com/view/?id=be9ac58a-0c5f-467e-8b5d-3cec677c3c95)と左下にあるExhaustivenessの値と思われるので、この値を100に合わせておきます。（直訳の通り、100回繰り返したのかもしれません。）

Exhaustivenessの値が高ければ高いほど、精度が高くなりますが、計算は重くなります。
デフォルト設定は8なので、最初に8~20あたりでスクリーニングすると良いでしょう。

![](https://storage.googleapis.com/zenn-user-upload/391d8021df7a-20241105.png)

右下のForwardをもしくは下のRunVinaを押すと分子ドッキングが始まります。

![](https://storage.googleapis.com/zenn-user-upload/ebf51c9c618d-20241105.png)

実行画面は以下のような感じです。

![](https://storage.googleapis.com/zenn-user-upload/6af95d550817-20241105.png)

ちなみにスクリーニング中の使用メモリは１GBくらいです。意外と少ないです。

![](https://storage.googleapis.com/zenn-user-upload/66975ad1039e-20241105.png)

一晩くらいかかるので、気軽に待ちましょう。

6000くらいでPyRxが落ちてしまうため、半分くらいで二回に分けるとよいでしょう。
