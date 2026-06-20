# DIA-MSプロテオミクス解析 よくある質問（FAQ）

## 🚀 環境構築・実行関連

### Q1. 「conda環境が見つからない」というエラーが出ます
**A.** 以下の順番で確認してください：

```bash
# 1. 環境が作成されているか確認
micromamba env list

# 2. crc-proteomics環境がない場合は作成
micromamba env create -f environment.yml

# 3. 環境をアクティベート
micromamba activate crc-proteomics

# 4. Pythonのバージョン確認
python --version  # 3.11.xが表示されるはず
```

### Q2. 「sage コマンドが見つかりません」エラーの対処法
**A.** Sageが正しくインストールされているか確認：

```bash
# 環境内でのインストール確認
micromamba activate crc-proteomics
sage --version

# インストールされていない場合
micromamba install -c bioconda sage-proteomics
```

### Q3. メモリ不足で処理が止まります
**A.** システムリソースを確認：

**推奨スペック**:
- **RAM**: 8GB以上（16GB推奨）
- **CPU**: 4コア以上
- **ストレージ**: 20GB以上の空き容量

**対処法**:
```python
# ファイル数を削減してテスト実行
mzml_files = mzml_files[:4]  # 最初の4ファイルのみ
```

### Q4. WSL環境での実行注意点
**A.** WSL特有の問題と対策：

```bash
# メモリ制限の確認
free -h

# 必要に応じて.wslconfigでメモリを増加
# C:\Users\[ユーザー名]\.wslconfig に追記:
# [wsl2]
# memory=8GB
```

---

## 📊 データ処理関連

### Q5. ダウンロードが途中で止まってしまいます
**A.** ネットワーク環境を確認：

**企業ネットワークの場合**:
```python
import requests
proxies = {
    'http': 'http://proxy.company.com:8080',
    'https': 'http://proxy.company.com:8080'
}
requests.get(url, proxies=proxies)
```

**回線が不安定な場合**:
- WiFiから有線LANに切り替え
- 混雑していない時間帯に実行
- 分割ダウンロードを検討

### Q6. Log2変換済みかどうかの判定方法
**A.** データの特徴で判定：

```python
def is_log_transformed(df):
    median_val = df.median().median()
    max_val = df.max().max()

    print(f"中央値: {median_val:.2f}")
    print(f"最大値: {max_val:.2f}")

    # 判定基準
    if median_val < 50 and max_val < 100:
        return True  # Log変換済み
    else:
        return False  # 生データ
```

### Q7. 欠損値補完のパラメータは変更すべきですか？
**A.** 実験条件に応じて調整：

**推奨設定**:
```python
# 標準設定（Toyota論文準拠）
IMPUTATION_DOWNSHIFT = 2.4
IMPUTATION_WIDTH = 0.3

# より保守的（信頼性重視）
IMPUTATION_DOWNSHIFT = 3.0
IMPUTATION_WIDTH = 0.2

# より積極的（検出力重視）
IMPUTATION_DOWNSHIFT = 2.0
IMPUTATION_WIDTH = 0.4
```

### Q8. 「FileNotFoundError」でスクリプトが止まります
**A.** ファイルパスを確認：

```bash
# 現在の作業ディレクトリを確認
pwd

# 期待されるファイル構造を確認
ls -la data/raw/
ls -la results/

# ファイルが存在するかチェック
find . -name "*.mzML" -type f
find . -name "human_proteome.fasta" -type f
```

---

## 🔬 プロテオミクス解析関連

### Q9. Sage vs DIA-NN、どちらを使うべき？
**A.** 用途に応じて選択：

| 項目 | Sage | DIA-NN |
|------|------|--------|
| **商用利用** | ✅ 制限なし | ❌ 制限あり |
| **検出感度** | 高い | より高い |
| **処理速度** | 高速 | より高速 |
| **安定性** | 安定 | やや不安定 |
| **推奨用途** | 企業研究・受託 | アカデミア |

### Q10. 検出タンパク質数が論文より少ないのはなぜ？
**A.** 複数の要因が考えられます：

**主な原因**:
1. **手法の違い**: Sage vs DIA-NN
2. **データ品質**: mzML変換での情報ロス
3. **パラメータ**: 検索設定の差
4. **サンプル数**: 論文では一部サンプルのみ使用の可能性

**改善策**:
```python
# 検索パラメータを緩める
"precursor_tol": "15 ppm",  # 10 -> 15に拡大
"fragment_tol": "15 ppm",   # 10 -> 15に拡大
```

### Q11. p値 < 0.05なのにqvalue > 0.05になるのはなぜ？
**A.** 多重検定補正の効果です：

```python
# 例：10,000回検定を実行
p_values = [0.01, 0.02, 0.03, ...]  # 生のp値
q_values = multipletests(p_values, method='fdr_bh')[1]  # 補正後

# p値0.03 → qvalue 0.08のような変化が起こる
# これは正常で、偽陽性率をコントロールするための調整
```

### Q12. PCAの寄与率が低いときはどう解釈すればよい？
**A.** データの性質を理解：

**寄与率の目安**:
- **PC1 > 30%**: 明確な主要パターン存在
- **PC1: 10-30%**: 複数要因が混在
- **PC1 < 10%**: 非常に複雑な構造

**対処法**:
```python
# より多くの主成分を確認
pca = PCA(n_components=10)
explained_ratio = pca.explained_variance_ratio_
cumsum_ratio = np.cumsum(explained_ratio)

# 累積寄与率で全体像を把握
print(f"PC1-3累積: {cumsum_ratio[2]:.2%}")
print(f"PC1-5累積: {cumsum_ratio[4]:.2%}")
```

---

## 📈 統計解析関連

### Q13. t検定でNaNが出る原因は？
**A.** データの問題を確認：

**よくある原因**:
```python
# 1. 分散がゼロ（すべて同じ値）
group1 = [5.0, 5.0, 5.0, 5.0]  # 分散=0
group2 = [6.0, 6.1, 5.9, 6.0]
# → t統計量が計算不可

# 2. サンプル数不足
group1 = [5.0]  # n=1
group2 = [6.0]  # n=1
# → 統計検定が不可能

# 3. 欠損値の残存
group1 = [5.0, np.nan, 6.0]  # NaNが残っている
```

**対処法**:
```python
# 事前チェックを追加
def safe_ttest(group1, group2):
    # NaN除去
    g1 = group1.dropna()
    g2 = group2.dropna()

    # サンプル数チェック
    if len(g1) < 2 or len(g2) < 2:
        return np.nan, np.nan

    # 分散チェック
    if g1.var() == 0 and g2.var() == 0:
        return np.nan, np.nan

    return stats.ttest_ind(g1, g2)
```

### Q14. 多重検定補正はなぜ必要？
**A.** 偽陽性の蓄積を防ぐため：

**例**: 10,000回の検定を実行
```python
# α=0.05で10,000回検定
# 偽陽性の期待値 = 10,000 × 0.05 = 500個

# 補正なし：500個の偽陽性が混入
# Benjamini-Hochberg法：50個以下に抑制（FDR=5%の場合）
```

### Q15. Volcano plotで有意なタンパク質が見つからない
**A.** 閾値設定を確認：

```python
# 厳しすぎる閾値の例
significant = (qvalues < 0.001) & (abs(log2fc) > 2.0)

# より現実的な閾値
significant = (qvalues < 0.05) & (abs(log2fc) > 1.0)

# さらに緩い探索的設定
significant = (pvalues < 0.05) & (abs(log2fc) > 0.5)
```

---

## 🛠️ トラブルシューティング

### Q16. スクリプト実行時のパーミッションエラー
**A.** ファイル権限を確認：

```bash
# 実行権限を付与
chmod +x scripts/run_all.sh

# スクリプトを実行
bash scripts/run_all.sh

# またはPythonで直接実行
python scripts/step_02_download_data.py
```

### Q17. Jupyter Notebookでカーネルが死ぬ
**A.** メモリとカーネル設定を確認：

```bash
# メモリ使用量確認
htop

# Jupyterを再起動
jupyter notebook stop 8888
jupyter notebook --port=8889 --no-browser
```

**Notebookの設定**:
```python
# セルの先頭に追加
%load_ext memory_profiler
%memit

# データサイズを制限
df_subset = df.sample(n=1000)  # 1000サンプルに削減
```

### Q18. 図の日本語が文字化けする
**A.** フォント設定を確認：

```python
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'DejaVu Sans'  # 英語フォント
# または
plt.rcParams['font.family'] = 'Noto Sans CJK JP'  # 日本語対応
```

---

## 📚 学習・理解関連

### Q19. プロテオミクス初心者におすすめの学習リソース
**A.** 段階的に学習：

**Step 1: 基礎概念**
- AJACS LC-MS解説動画（前編・後編）
- 本プロジェクトの [glossary.md](glossary.md)

**Step 2: 実践**
- 本プロジェクトの各スクリプト解説
- 少数ファイルでの試行

**Step 3: 発展**
- Nature Protocols誌のプロテオミクス論文
- Bioconductorのプロテオミクスパッケージ

### Q20. 結果をどう論文にまとめればよい？
**A.** 標準的な構成に従って：

**Methods**:
```
DIA-MS analysis was performed using Sage-proteomics v0.14
with the following parameters: enzyme = trypsin,
missed cleavages ≤ 1, peptide length 7-45 AA,
precursor tolerance = 10 ppm, fragment tolerance = 10 ppm,
FDR < 1%.
```

**Results**:
- 検出タンパク質数
- 有意差タンパク質数
- 主要なバイオマーカー候補

**Figures**:
- PCAプロット（サンプルのクラスタリング）
- Volcanoプロット（有意差の可視化）
- ヒートマップ（主要タンパク質の発現パターン）

---

💡 **さらに詳しい情報が必要な場合**:
1. 各スクリプトの詳細解説（`explanations/step_XX_*.md`）
2. 用語集（`explanations/glossary.md`）
3. プロジェクト全体概要（`explanations/00_overview.md`）

を参照してください。