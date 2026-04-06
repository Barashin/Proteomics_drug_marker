---
title: "COSMICデータベースでがん関連タンパク質を特定する【論文再現シリーズ #8】"
emoji: "🎯"
type: "tech"
topics: ["proteomics", "cancer", "COSMIC", "labcode"]
published: false
---

# COSMICデータベースでがん関連タンパク質を特定する

## はじめに

DIA-MS で同定したタンパク質の中に、がんに関連するタンパク質はどれくらい含まれているか？論文では10,329 タンパク質のうち 531 個（71%）が COSMIC がん関連として特定されました。本書のサブセット（sage + 18ファイル、2,234 タンパク質）では **21 個（全がん関連 10.5%）、CRC特異的 7 個（10.8%）** を検出。この記事では、**COSMIC（Catalogue Of Somatic Mutations In Cancer）** データベースとの照合手順を解説します。

## 前提

- [#7 差分発現解析](article-07-differential.md) が完了していること

## スクリプト全文: step_05_cosmic_analysis.py

> コアとなる照合ロジックを抜粋します。完全版は `scripts/step_05_cosmic_analysis.py` を参照。

```python
#!/usr/bin/env python3
"""
Step 5: COSMICデータベースとの照合
===================================
同定タンパク質をCOSMIC Cancer Gene Censusと照合し、
がん関連タンパク質のカバー率を計算するスクリプトです。

【COSMICとは？】
  Catalogue Of Somatic Mutations In Cancer の略。
  がんの体細胞変異に関する世界最大のデータベース。
  Cancer Gene Census: がんに繰り返し変異が見つかる遺伝子（タンパク質）のリスト。
  - 全がん種: 748遺伝子
  - 大腸がん特異的: 65遺伝子
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# COSMICのがん関連遺伝子リスト
# ============================================================
# Cancer Gene Census から取得した代表的ながん関連遺伝子
# 本来は748遺伝子だが、ここでは代表的なものを掲載

# 大腸がん（CRC）に特異的なCOSMIC遺伝子（65遺伝子のうち代表的なもの）
COSMIC_CRC_GENES = [
    "APC",     # Wntシグナル経路の腫瘍抑制因子。CRCの約80%で変異
    "TP53",    # がん抑制遺伝子の代表格。細胞周期停止・アポトーシスを制御
    "KRAS",    # RAS-MAPK経路のがん遺伝子。CRCの約40%で変異
    "BRAF",    # MAPK経路のキナーゼ。V600E変異がCRCの10%に
    "PIK3CA",  # PI3K経路の触媒サブユニット。増殖シグナルに関与
    "SMAD4",   # TGF-βシグナル経路。転移に関連
    "FBXW7",   # ユビキチンリガーゼ。Notch, MYC等の分解を制御
    # ... （省略: 完全版はスクリプト内に記載）
]

# 全がん種のCOSMIC遺伝子（代表的なサブセット）
COSMIC_ALL_GENES = COSMIC_CRC_GENES + [
    "ABL1", "AKT1", "ALK", "AR", "ARID1A", "ATM", "BCL2",
    "BRCA1", "BRCA2", "CDK4", "CDKN2A", "EGFR", "ERBB2",
    "EZH2", "FGFR1", "FLT3", "IDH1", "JAK2", "KIT", "MET",
    "MYC", "NF1", "NOTCH1", "PTEN", "RAF1", "RB1", "RET",
    "STAT3", "VHL", "WT1",
    # ... （省略: 完全版はスクリプト内に記載）
]


# ============================================================
# COSMIC照合のコアロジック
# ============================================================
def cosmic_analysis(identified_proteins):
    """
    同定タンパク質とCOSMIC Cancer Gene Censusを照合する関数。

    【処理内容】
      1. 同定タンパク質のリストとCOSMICリストの共通要素（積集合）を計算
      2. カバー率 = 共通要素数 / COSMIC登録数 × 100 を算出

    引数:
        identified_proteins (set): 同定タンパク質名の集合

    戻り値:
        dict: 照合結果（カバー率、共通タンパク質リスト等）
    """
    cosmic_all = set(COSMIC_ALL_GENES)     # 全がん関連（748遺伝子の代表）
    cosmic_crc = set(COSMIC_CRC_GENES)     # CRC特異的（65遺伝子の代表）

    # 積集合（&演算子）: 両方に含まれるタンパク質を取得
    overlap_all = identified_proteins & cosmic_all
    overlap_crc = identified_proteins & cosmic_crc

    return {
        "cosmic_all_total": len(cosmic_all),
        "cosmic_crc_total": len(cosmic_crc),
        "overlap_all": overlap_all,
        "overlap_crc": overlap_crc,
        # カバー率の計算
        "coverage_all": len(overlap_all) / len(cosmic_all) * 100,
        "coverage_crc": len(overlap_crc) / len(cosmic_crc) * 100,
    }
```

---

## コード詳細

### 集合演算（set）の使い方

```python
# 積集合: 両方に含まれる要素
overlap = set_A & set_B

# 差集合: Aにはあるが、Bにはない要素
only_in_A = set_A - set_B
```

- Pythonの `set` 型は集合演算が得意で、数万件のリスト照合も一瞬で完了します
- リスト（`list`）で同じことをすると `O(n²)` かかりますが、`set` なら `O(n)` です

### 論文と本書の結果比較

| 項目 | 論文 (DIA-NN, 10329タンパク質) | 本書 (sage, 2234タンパク質) |
|------|-----------------------------|-----------------------------|
| 全がん関連カバー率 | 531/748 = **71%** | 21/200 = **10.5%** |
| CRC特異的カバー率 | 48/64 = **75%** | 7/65 = **10.8%** |
| 検出された主要ドライバー | - | **KRAS, CTNNB1, PIK3CA 等 21個** |

### COSMICカバレッジの可視化

![COSMICカバレッジ](images/fig_cosmic_coverage.png)

COSMIC Cancer Gene Censusに登録されているがん関連遺伝子のうち、本書のDIA-MSプロテオミクスで同定できた割合を棒グラフで示しています。全がん関連では200個中21個（10.5%）、CRC特異的では65個中7個（10.8%）をカバーしています。

### カバー率が論文より低い理由

| 理由 | 説明 |
|------|------|
| 同定タンパク質数の差 | sage + 18ファイル = 2,234 vs DIA-NN + 32ファイル = 10,329 |
| 低発現分子の検出感度 | DIA-NN の深層学習ライブラリ予測が低発現タンパク質に強い |
| MBR 未実装 | sage 0.14 は Match Between Runs 未対応 |

ただし、**KRAS などの大腸がん主要ドライバー遺伝子は本書でも検出** できています。カバー率の絶対値は小さくても、バイオロジカルに重要な遺伝子は捕捉できているということです。

## 実行方法

```bash
python scripts/step_07_cosmic_analysis.py
```

## まとめ

DIA-MSプロテオミクスは、がん関連タンパク質をカバーする強力な手法です。論文では DIA-NN で約 70-75% をカバーし、本書の sage + 縮小サブセットでは約 10-11% に留まりますが、**主要なドライバー遺伝子（KRAS, CTNNB1, PIK3CA 等）は商用クリアな sage パイプラインでも検出** できることが確認できました。

> 前回: [#7 差分発現解析](article-07-differential.md)
> 次回: [#9 ステージ別解析](article-09-stage.md) — ANOVAとクラスター分析

#バイオインフォマティクス #プロテオミクス #COSMIC #がん研究 #labcode
