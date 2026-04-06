---
title: "📘 第12章：RMSD, RMSF, RoG, hydrogen bond, RDF解析"
free: false
---

本章はリガンド-タンパク質 MDシミュレーション（分子動力学シミュレーション）を論文に沿って行なった結果に対し、
論文で使われている解析方法であるRMSD, RMSF, RoG, hydrogen bond, RDF解析を行なっています。

以下のような図がRMSD, RMSF, RoG, hydrogen bond, RDF解析の結果作成できます。
![](https://storage.googleapis.com/zenn-user-upload/6d1263c74187-20241125.png)

```
動作検証済み環境

Windows 11 Home, 13th Gen Intel(R) Core(TM) i7-13700, 
64 ビット オペレーティング システム、x64 ベース プロセッサ, メモリ：32GB
```

本章では論文のFigure5, 7の作り方を説明します。
分子ドッキングの結果が違うので、正確に同じグラフにはなりませんが、作成の仕方を参考にしてみてください。

## 参考文献
https://www.udemy.com/share/105vCm3@9eQi5TNCBbh2arjqlQq11eIj5xU4FyZWgdCUVdVrt0drqyT67M0pVeRk39zckzTBVw==/