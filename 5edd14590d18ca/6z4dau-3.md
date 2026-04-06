---
title: "📰 結果"
free: false
---

**PyRxは結構落ちやすいソフトなので、結果が出たら、早めにsaveしてください。結果は最初に設定したOutputフォルダにあるので、途中で落ちた場合は途中のリガンドから再開すると良いでしょう。**

結果の右上のsaveボタンでsaveできます。

Binding Energyの所を箇所を押して、エネルギーが低い順に並び替えます。

![](https://storage.googleapis.com/zenn-user-upload/3c0b5dd396a2-20241105.jpg)

E（エネルギー）が高すぎるものについては、エネルギー最小化のところでバグっているので、取り除きます。

CMNPD14737のもの

実際開いてみるとバグっていることがわかります。

![](https://storage.googleapis.com/zenn-user-upload/042061b3761e-20241105.png)

エネルギー計算がうまくいっていないものを取り除くと、

最もbinding affinityが高いのはCMNPD14239で、-7.7 kcal/molになりました。

![](https://storage.googleapis.com/zenn-user-upload/db01ae1a2481-20241105.jpg)

ちなみに、論文で紹介されているCMNPD28986については-5.5kcal/mol(論文では-7.9 kcal/mol)になりました。

選択して、Moleculeタブに行くと、画面に結合を移すことができます。

![](https://storage.googleapis.com/zenn-user-upload/0cdf0f9506e8-20241105.jpg)

詳細な解析については他の可視化ツールの方がよいので、そちらで可視化させます。

一番結合力の高いリガンドのPDBファイルを保存してください。

![](https://storage.googleapis.com/zenn-user-upload/9255928cb269-20241105.jpg)

この保存したものとタンパク質ファイルを第1章で紹介したPyMOLなどで開いてみてください。

タンパク質-リガンドの可視化については、
**第8章：タンパク質-リガンドの相互作用可視化**
で用いるPLIP, PoseViewの方が綺麗に見れます。こちらをご参照ください。
