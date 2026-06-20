---
name: article-editor
description: 高度記事エディタ - ワークフロー最終品質管理ステップ。book-designで生成された技術書/ブログ全体の構成・整合性・図表の包括的チェックを行い、問題があればbook-designに具体的フィードバックを返すゲートキーパー機能。記事全体の論理的流れ、章間の矛盾・重複、図表の視覚的整合性、用語統一、参照整合性を多角的に分析し、出版品質まで引き上げる。「記事の最終チェック」「全体構成チェック」「記事エディット」「品質チェック」「整合性チェック」「article-editor」「final review」「構成見直し」等でトリガーする。
---

# Article Editor - 高度記事エディタ

記事全体の品質を多角的に分析し、出版レベルまで引き上げるワークフロー最終ゲートキーパー。

## 機能概要

### 1. 記事構成の包括的分析
- **構造整合性**: 記事番号連続性、章立て一貫性
- **内容整合性**: タイトル重複、論理的流れ、セクション構造
- **参照整合性**: 内部リンク、相互参照の正確性
- **記事バランス**: 記事長の適切性、内容分散
- **記事統合推奨**: 細分化記事の統合候補検出、Zenn出版最適化

### 2. 技術的整合性チェック
- **用語統一**: 専門用語の表記統一性
- **図表整合性**: 図番号、キャプション、参照の一貫性
- **コード一貫性**: サンプルコード、パス参照の整合性

### 3. 図の視覚的品質管理
- **存在確認**: 画像ファイルの実在性
- **品質評価**: Alt text、キャプションの充実度
- **論文比較**: 元論文図との視覚的整合性（可能な場合）

### 4. 修正フィードバック機能
- **問題の分類**: エラー・警告・情報レベルでの問題分類
- **book-design連携**: 具体的修正指示でbook-designに差し戻し
- **改善提案**: 構成改善・統合・分割の提案

## 実行手順

### Step 1: 記事全体の構成分析

スクリプト `scripts/check_article_structure.py` を実行して記事構成を分析：

```bash
python scripts/check_article_structure.py Article/blog/
```

**チェック項目：**
- ファイル命名規則の一貫性
- 記事番号の連続性と論理性
- タイトル重複・類似性
- セクション構造パターン
- 内部リンクの有効性
- 記事長のバランス

### Step 2: 図表の整合性分析

スクリプト `scripts/check_figure_consistency.py` で図表を分析：

```bash
# 基本的な図チェック
python scripts/check_figure_consistency.py Article/blog/

# 論文図と比較する場合
python scripts/check_figure_consistency.py Article/blog/ data/paper_figures/
```

**チェック項目：**
- 図ファイルの存在確認
- 図番号と実際の図数の整合性
- Alt textの品質評価
- 重複画像の検出
- 論文図との類似性比較

### Step 3: 総合品質評価

両スクリプトの結果を統合して総合評価を実施：

**エラーレベル（Critical）:**
- 存在しない画像への参照
- 重複するタイトル
- 破綻した内部リンク
→ **book-designに即座に差し戻し**

**警告レベル（Warning）:**
- 記事番号の飛び
- 用語の表記ゆれ
- 細分化記事の統合推奨（02a,02b等）
- 不十分なAlt text
→ **改善推奨、必要に応じて修正**

**情報レベル（Info）:**
- 類似構成記事の統合提案
- 論文図との比較推奨
→ **最適化の提案**

### Step 4: フィードバック生成

問題に応じて適切なアクションを決定：

#### 4.1 book-designへの差し戻し基準

以下の場合はbook-designに差し戻し：
- **構造的問題**: 5個以上のエラーレベル問題
- **整合性破綻**: 重複タイトル、破綻リンク
- **論理的矛盾**: 章の順序が不適切
- **技術的エラー**: 存在しない図への参照

#### 4.2 差し戻し時のフィードバック内容

```markdown
## book-design修正指示

### 重要な問題（即修正要）
- [具体的エラー内容]
- 影響ファイル: [ファイル名]
- 推奨修正: [具体的な修正方法]

### 構成改善提案
- [記事統合/分割提案]
- [章立て見直し提案]
- [用語統一提案]

### 図表関連
- [図の修正/追加提案]
- [キャプション改善]
```

#### 4.3 軽微な修正（自己完結）

警告・情報レベルは記事エディタ内で修正：
- Alt textの改善
- 用語統一の軽微な修正
- 図番号の調整

### Step 5: 品質基準の判定

**出版可能基準：**
- エラーレベル問題：0件
- 警告レベル問題：2件以下
- 記事構成の論理的一貫性確保
- 図表の完全性確保

基準を満たさない場合は該当ステップに差し戻し。

## 記事構成最適化（新機能）

### 記事統合推奨機能

**方針**:
- 文字数制限は緩和（目的・機能優先）
- 細分化記事（02a,02b等）の統合を推奨
- Zenn出版時のチャプター構成を考慮

**統合候補の自動検出**:
```python
def detect_consolidation_candidates(articles):
    """統合推奨記事を検出"""

    # 同一番号の細分化記事
    numbered_groups = group_by_number(articles)

    suggestions = []
    for number, group in numbered_groups.items():
        if len(group) > 1:  # 02a, 02b等の細分化
            # 機能名を推定（最初の記事から）
            function_name = extract_function_name(group[0])
            target_article = f"article-{number:02d}-{function_name}.md"
            target_notebook = f"notebook-{number:02d}-{function_name}.ipynb"

            suggestions.append({
                "type": "consolidation",
                "articles": group,
                "target_article": target_article,
                "target_notebook": target_notebook,
                "reason": "機能単位での統合推奨（記事とnotebook名完全対応）"
            })

    return suggestions
```

**Zennチャプター構成提案**:
- 統合記事 → 自然な「## セクション」構成
- 各セクション → Zenn出版時の小見出し
- 読みやすい流れ → 学習効率向上

## 使用例

### 基本的な品質チェック

```bash
# 記事全体の構成チェック
python scripts/check_article_structure.py Article/blog/ > structure_report.json

# 図の整合性チェック
python scripts/check_figure_consistency.py Article/blog/ > figure_report.json
```

### 問題に応じた対処

**Case 1: 構成問題が多数**
→ book-designに差し戻し、記事構成の根本的見直し

**Case 2: 図の問題のみ**
→ 図の修正・Alt text改善を実施

**Case 3: 軽微な修正のみ**
→ 記事エディタ内で修正完了

## 高度な機能

### カスタム検証ルール

プロジェクト固有の検証ルールを `references/custom_rules.md` で定義可能：

- 特定用語の強制統一
- 図の命名規則
- セクション構造の規約

### 継続的品質管理

定期的な品質チェックで品質維持：

```bash
# 週次品質レポート生成
python scripts/weekly_quality_check.py Article/blog/
```

### 他ステップとの連携

- **code-debug**: 図生成コードの修正指示
- **book-design**: 構成・内容修正指示
- **publication-editor**: 最終出版判定への品質保証

## トラブルシューティング

### よくある問題

**記事数が多すぎて管理困難**
→ 記事統合の提案生成

**図の品質が一定でない**
→ 図生成コードの標準化提案

**用語の表記ゆれが多い**
→ 用語集の作成・統一処理の自動化

記事エディタは品質のゲートキーパーとして、完璧な記事品質を保証し、読者に価値ある技術コンテンツを提供する最後の砦となる。