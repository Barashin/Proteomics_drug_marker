---
name: zenn-publisher
description: Zenn出版最適化スキル - 統合記事をZenn bookとして出版するための最終調整を行う。記事統合、チャプター構成最適化、小見出し調整、技術書体裁の整備、Zenn固有の機能（コードブロック、注釈、リンク）への対応。labcode-publisherワークフローの最終ステップとして、読みやすく完成度の高いZenn技術書を生成する。「Zenn出版」「チャプター調整」「技術書最適化」「zenn-publisher」等でトリガー。
---

# Zenn Publisher - Zenn出版最適化

記事をZenn book形式に最適化し、読みやすい技術書として出版する。

## 機能概要

### 1. 記事統合とチャプター構成
- **統合記事対応**: 細分化記事を機能単位でまとめ
- **チャプター最適化**: Zenn bookの理想的な章立て
- **小見出し調整**: 読みやすい見出しレベルの設定

### 2. Zenn固有機能の活用
- **コードブロック最適化**: 言語指定、ファイル名表示
- **注釈ブロック活用**: :::message, :::details, :::alert
- **図表の最適化**: Zennでの表示に適した形式

### 3. 技術書体裁の整備
- **目次生成**: 自動的な章立て構成
- **相互リンク**: 章間の参照関係
- **読者ガイド**: 前提知識、環境要件の明確化

## 実行手順

### Step 1: 記事統合分析

統合対象記事の分析と最適化方針の決定：

```python
def analyze_consolidation_targets(blog_dir):
    """統合対象記事を分析"""

    articles = scan_articles(blog_dir)

    # 統合グループの特定
    consolidation_groups = {
        "data": ["02a-data-acquisition", "02b-data-formats"],
        "convert": ["03a-convert-basics", "03b-convert-execution"],
        "sage": ["04a-sage-fundamentals", "04b-sage-execution",
                "04c-sage-gene-mapping", "04d-sage-protein-aggregation"],
        # ...
    }

    # Zenn章構成の提案
    zenn_chapters = generate_zenn_structure(consolidation_groups)

    return {
        "consolidation_plan": consolidation_groups,
        "zenn_structure": zenn_chapters,
        "optimization_targets": identify_optimization_needs(articles)
    }
```

### Step 2: Zennチャプター構成最適化

理想的なZenn book構造の生成：

**推奨チャプター構成**:
```
Zenn Book: "DIA-MSプロテオミクスで大腸がんバイオマーカー探索"

├── Chapter 0: はじめに（概要・環境構築）
├── Chapter 1: データ取得と変換
│   ├── 1.1 ProteomeXchangeからのデータ取得
│   ├── 1.2 mzMLファイル形式の理解
│   └── 1.3 RAW→mzML変換
├── Chapter 2: Sage-proteomicsによるDIA解析
│   ├── 2.1 Sage基礎と環境設定
│   ├── 2.2 DIA検索実行
│   ├── 2.3 遺伝子マッピング
│   └── 2.4 タンパク質マトリクス構築
├── Chapter 3: データ探索と前処理
│   ├── 3.1 データ分布の確認
│   ├── 3.2 相関分析
│   └── 3.3 欠損値処理
├── Chapter 4: 可視化と解析
│   ├── 4.1 基本的な可視化
│   ├── 4.2 主成分分析
│   └── 4.3 差分発現解析
└── Chapter 5: 高度な解析（OpenMS編）
    ├── 5.1 OpenMS深層学習解析
    ├── 5.2 ステージ別解析
    └── 5.3 手法比較
```

### Step 3: 記事内容の統合とセクション調整

```python
def create_unified_chapter(source_articles, chapter_config):
    """複数記事を1つのZennチャプターに統合"""

    chapter_content = []

    # チャプター導入
    chapter_content.append(generate_chapter_intro(chapter_config))

    # 各セクションの統合
    for i, (article, section_config) in enumerate(zip(source_articles, chapter_config['sections'])):

        # セクション見出し（## レベルに統一）
        section_header = f"## {section_config['number']} {section_config['title']}"
        chapter_content.append(section_header)

        # セクション概要
        if section_config.get('description'):
            chapter_content.append(f"*{section_config['description']}*\n")

        # 記事内容の取り込み（見出しレベル調整）
        article_content = load_article_content(article)
        adjusted_content = adjust_heading_levels(article_content, base_level=3)
        chapter_content.append(adjusted_content)

        # セクション間の区切り
        if i < len(source_articles) - 1:
            chapter_content.append("\n---\n")

    # チャプターまとめ
    chapter_content.append(generate_chapter_summary(chapter_config))

    return "\n".join(chapter_content)
```

### Step 4: Zenn固有機能の最適化

**注釈ブロックの活用**:
```markdown
:::message
本シリーズに入る前に、LC-MSの基礎知識を確認することをお勧めします。
:::

:::details 上級者向け設定
より高度な設定については以下を参照してください。
:::

:::alert
この操作は大量のメモリを使用します。16GB以上のRAMを推奨します。
:::
```

**コードブロックの最適化**:
```python
def optimize_code_blocks(content):
    """Zenn向けコードブロック最適化"""

    # 言語指定の追加
    content = re.sub(r'```\n', '```python\n', content)

    # ファイル名の表示
    content = re.sub(
        r'```python\n(# === .+ ===.*?)\n',
        r'```python:example.py\n\1\n',
        content,
        flags=re.DOTALL
    )

    # 実行結果の明示
    content = add_execution_results_blocks(content)

    return content
```

### Step 5: 技術書品質の向上

**相互参照の最適化**:
```markdown
[第2章で説明したSage設定](./02-sage-analysis)を使用して、実際のデータ解析を行います。

実行結果の詳細な解釈については、[第4章の可視化手法](./04-visualization)も併せてご覧ください。
```

**読者ガイドの充実**:
```markdown
## この章で学ぶこと
- [ ] ProteomeXchangeからのデータ取得方法
- [ ] mzMLファイル形式の理解
- [ ] RAW→mzML変換の実践

## 前提条件
- Python 3.11以上
- 16GB以上のRAM
- 第1章の環境構築完了

## 実行時間目安
約30分（データダウンロード時間を除く）
```

## 高度な機能

### 目次とナビゲーションの自動生成

```python
def generate_book_navigation(chapters):
    """Zenn book全体のナビゲーション生成"""

    # config.yaml生成
    config = {
        "title": "DIA-MSプロテオミクスで大腸がんバイオマーカー探索",
        "summary": "無料ツールのみでDIA-MSプロテオミクス解析を完全再現",
        "topics": ["proteomics", "bioinformatics", "python", "labcode"],
        "published": False,
        "price": 0  # 無料公開
    }

    # chapters.yaml生成
    chapters_config = []
    for chapter in chapters:
        chapters_config.append({
            "title": chapter["title"],
            "file": chapter["filename"]
        })

    return config, chapters_config
```

### Zenn固有の最適化

**画像の最適化**:
```python
def optimize_images_for_zenn(content, images_dir):
    """Zenn向け画像最適化"""

    # 画像サイズの調整
    optimized_images = resize_images_for_web(images_dir)

    # 画像パスの修正
    for img_path in optimized_images:
        rel_path = f"/images/{img_path.name}"
        content = content.replace(f"images/{img_path.name}", rel_path)

    return content
```

**外部リンクの整理**:
```python
def organize_external_links(content):
    """外部リンクの整理と注釈"""

    # 重要なリンクは注釈ブロックで強調
    content = re.sub(
        r'\[([^\]]+)\]\((https://github\.com/[^)]+)\)',
        r':::message\n[\1](\2)\n:::\n',
        content
    )

    return content
```

## 品質チェックとZenn準拠

### Zenn品質基準

**必須要件**:
- [ ] 各章に明確な学習目標
- [ ] コードブロックに適切な言語指定
- [ ] 図表に代替テキスト
- [ ] 外部依存の最小化

**推奨要件**:
- [ ] 注釈ブロックの効果的活用
- [ ] 章間の自然な流れ
- [ ] 実行可能なコード例
- [ ] トラブルシューティング情報

### 自動品質チェック

```python
def validate_zenn_compliance(book_content):
    """Zenn準拠性の自動チェック"""

    issues = []

    # コードブロックの言語指定チェック
    unspecified_blocks = re.findall(r'```\n(?!.*```)', book_content)
    if unspecified_blocks:
        issues.append("言語指定なしのコードブロックがあります")

    # 章立ての一貫性チェック
    headings = re.findall(r'^#+\s+(.+)$', book_content, re.MULTILINE)
    if not validate_heading_hierarchy(headings):
        issues.append("見出し階層が不正です")

    return issues
```

## 使用例

### 基本的なZenn book生成

```bash
# 記事統合からZenn book生成まで
python scripts/generate_zenn_book.py \
  --articles Article/blog/ \
  --output zenn-book/ \
  --config zenn_config.yaml
```

### 個別チャプター最適化

```bash
# 特定章のみZenn最適化
python scripts/optimize_for_zenn.py \
  --input article-04-sage.md \
  --output zenn-book/04-sage-analysis.md \
  --chapter-config chapter_04.yaml
```

## 他スキルとの連携

### labcode-publisherとの統合

**Step A8.5: zenn-publisher**（article-editor後に挿入）
- article-editorで品質確認済み記事をZenn最適化
- book-designの成果物をZenn book形式に変換
- publication-editorでZenn品質基準での最終判定

### 統合ワークフロー例

```
book-design → article-editor → zenn-publisher → publication-editor
     ↓              ↓              ↓                  ↓
   記事作成      品質チェック    Zenn最適化        最終出版判定
```

Zenn Publisherは技術記事を読者フレンドリーな技術書に変換し、Zenn固有の機能を最大限活用した高品質な出版物を生成する。