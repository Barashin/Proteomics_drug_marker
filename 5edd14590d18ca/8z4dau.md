---
title: "📘 第11章：リガンド-タンパク質 MDシミュレーション（分子動力学シミュレーション）"
free: false
---

本章はリガンド-タンパク質 MDシミュレーション（分子動力学シミュレーション）を論文に沿って、行っています。分子のエネルギーや運動を計算するための相互作用モデルである力場はff14SBを使用しており、[GROMACSのチュートリアル](https://labo-code.com/bioinformatics/gromax-protein-small-molecule-first-part/)とは異なる方法で行なっています。
エラーの回避方法も含め、詳し目に書いているので、ぜひご覧ください！

```
動作検証済み環境

Windows 11 Home, 13th Gen Intel(R) Core(TM) i7-13700,
64 ビット オペレーティング システム、x64 ベース プロセッサ, メモリ：32GB
```

今回はMaterial and Methodの「EPSP-inhibitors dynamics analysis」の部分をフォローしています。

各パラメータは以下になります。

- 力場：AMBER20、Ff14SB(タンパク質）、GAFF（リガンド)
- 水モデル:TIP3P
- エネルギー最小化:1000ステップ
- 温度:徐々に310K、Langevinダイナミクスで20psの温度上昇
- NPT、NVT平衡化:100ps
- MD:100ns

論文ではAMBERで行っていますが、今回はGROMACSで行っています。適宜パラメータを変えています。

本章ができると、以下のようにsimulationできます。
![](https://storage.googleapis.com/zenn-user-upload/d45858872f7f-20241125.gif)

## 全体の流れ

本simulationは以下の流れで行っていきます。

1. リガンドの準備
2. タンパク質の準備
3. リガンド-タンパク質の準備
4. 水とイオン追加
5. エネルギー最小化
6. NVT（定温）とNPT（定圧定温）の平衡化
7. MDシミュレーション
8. 可視化


## 本章で用いられるツール一覧

- [GROMACS](https://www.gromacs.org/): MDシミュレーション各種設定
- [acpype](https://github.com/alanwilter/acpype): リガンドへの力場付与


## 環境構築

**第9章：MD simulationのための環境構築** で構築した環境をそのまま使用します。
まだ環境構築がお済みでない方は、先にそちらをご覧ください。

```bash
# 環境をアクティベート
micromamba activate gromacs_env
```


## 参考文献

https://www.udemy.com/share/105vCm3@9eQi5TNCBbh2arjqlQq11eIj5xU4FyZWgdCUVdVrt0drqyT67M0pVeRk39zckzTBVw==/
https://qiita.com/Ag_smith/items/6fd0afa8effc119b0f38
