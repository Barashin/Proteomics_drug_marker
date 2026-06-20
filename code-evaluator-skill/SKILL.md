---
name: code-evaluator
description: コード評価スキル - code-implementで実装したコードを実行前に静的・動的に評価し、複数の評価指標（構文・依存関係・インポート・パス参照・スタイル・再現性・スモーク実行）でエラーなく実行できる状態かを判定する。問題があれば具体的な修正提案付きで評価レポートを出力する。「評価」「コード評価」「事前チェック」「実行前チェック」「静的解析」「evaluate」「コードは動くか」「エラーはないか」「品質チェック」等でトリガーする。
---

# Code Evaluator Skill

code-implementで実装したJupyter Notebookコードを実行前に多角度から評価し、実行成功率を最大化する。

## 評価フレームワーク

### 7軸品質評価システム

| 評価軸 | 重要度 | 評価内容 | 合格基準 |
|-------|--------|----------|----------|
| **構文** | 🔴 Critical | Python構文エラー・インデントエラー | 0エラー |
| **依存関係** | 🔴 Critical | import文・ライブラリ互換性 | 全依存関係解決済み |
| **パス参照** | 🟡 Major | ファイルパス・ディレクトリ参照 | 相対パス・存在確認済み |
| **スタイル** | 🟢 Minor | PEP8準拠・可読性 | 90%準拠 |
| **再現性** | 🔴 Critical | 実行順序・状態依存性 | 上から順次実行可能 |
| **スモーク実行** | 🔴 Critical | 基本機能の動作確認 | エラーなし |
| **環境互換性** | 🟡 Major | OS・Python版・パッケージ版 | 主要環境で動作 |

## 評価ワークフロー

### Phase 1: 静的コード解析

**1.1 構文チェック**:
```python
def syntax_check(notebook_path):
    """Jupyter Notebookの全コードセルを構文チェック"""

    import ast
    import nbformat

    nb = nbformat.read(notebook_path, as_version=4)
    syntax_errors = []

    for cell_idx, cell in enumerate(nb.cells):
        if cell.cell_type == 'code':
            try:
                ast.parse(cell.source)
                print(f"✅ Cell {cell_idx}: Syntax OK")
            except SyntaxError as e:
                error_detail = {
                    'cell_index': cell_idx,
                    'line_number': e.lineno,
                    'error_message': str(e),
                    'code_snippet': cell.source.split('\n')[e.lineno-1] if e.lineno else ''
                }
                syntax_errors.append(error_detail)
                print(f"❌ Cell {cell_idx}: Syntax Error - {e}")

    return syntax_errors
```

**1.2 依存関係解析**:
```python
def dependency_analysis(notebook_path):
    """importされるパッケージと環境ファイルの整合性チェック"""

    import re
    from pathlib import Path

    # Notebookからimport文を抽出
    nb_imports = extract_imports_from_notebook(notebook_path)

    # environment.yml または requirements.txt を読み込み
    env_packages = load_environment_packages()

    missing_deps = []
    for imp in nb_imports:
        if imp not in env_packages:
            missing_deps.append(imp)

    # 各パッケージのバージョン互換性チェック
    version_conflicts = check_version_compatibility(env_packages)

    return {
        'missing_dependencies': missing_deps,
        'version_conflicts': version_conflicts,
        'total_imports': len(nb_imports),
        'resolved_imports': len(nb_imports) - len(missing_deps)
    }
```

**1.3 パス参照チェック**:
```python
def path_reference_check(notebook_path):
    """ファイルパス・ディレクトリ参照の妥当性チェック"""

    import re
    from pathlib import Path

    nb_content = read_notebook_as_text(notebook_path)

    # ファイルパスのパターンマッチング
    path_patterns = [
        r'["\']([^"\']+\.[a-zA-Z]+)["\']',  # "file.ext" or 'file.ext'
        r'Path\(["\']([^"\']+)["\']\)',     # Path("path")
        r'pd\.read_csv\(["\']([^"\']+)["\']',  # pd.read_csv("file")
        r'open\(["\']([^"\']+)["\']'        # open("file")
    ]

    referenced_paths = []
    for pattern in path_patterns:
        matches = re.findall(pattern, nb_content)
        referenced_paths.extend(matches)

    # パスの存在確認
    path_issues = []
    notebook_dir = Path(notebook_path).parent

    for path_str in referenced_paths:
        # 絶対パス・相対パス・プロジェクト内パスの判定
        path_obj = Path(path_str)

        if path_obj.is_absolute():
            if not path_obj.exists():
                path_issues.append({
                    'path': path_str,
                    'issue': 'absolute_path_not_found',
                    'suggestion': f'Check if {path_str} exists on target system'
                })
        else:
            # 相対パスの場合、notebook基準で存在確認
            full_path = notebook_dir / path_obj
            if not full_path.exists():
                path_issues.append({
                    'path': path_str,
                    'issue': 'relative_path_not_found',
                    'suggestion': f'Create {full_path} or adjust path reference'
                })

    return path_issues
```

### Phase 2: 動的コード評価

**2.1 実行順序依存性チェック**:
```python
def execution_order_check(notebook_path):
    """セル間の実行順序依存性を解析"""

    nb = nbformat.read(notebook_path, as_version=4)

    defined_variables = set()
    execution_issues = []

    for cell_idx, cell in enumerate(nb.cells):
        if cell.cell_type == 'code':
            # 変数定義の抽出
            new_vars = extract_variable_definitions(cell.source)

            # 変数使用の抽出
            used_vars = extract_variable_usage(cell.source)

            # 未定義変数の使用チェック
            undefined_vars = used_vars - defined_variables
            if undefined_vars:
                execution_issues.append({
                    'cell_index': cell_idx,
                    'issue': 'undefined_variables',
                    'variables': list(undefined_vars),
                    'suggestion': f'Define {undefined_vars} before this cell'
                })

            # 定義済み変数リストを更新
            defined_variables.update(new_vars)

    return execution_issues
```

**2.2 スモーク実行テスト**:
```python
def smoke_execution_test(notebook_path):
    """Notebookの基本実行テスト（実際には実行せず、実行可能性を評価）"""

    import tempfile
    import subprocess
    from pathlib import Path

    # テンプレートNotebookの作成（重要セルのみ抜粋）
    smoke_nb = create_smoke_test_notebook(notebook_path)

    with tempfile.NamedTemporaryFile(suffix='.ipynb', delete=False) as tmp_file:
        nbformat.write(smoke_nb, tmp_file.name)

        # jupyter nbconvert --execute でドライラン
        cmd = [
            'jupyter', 'nbconvert',
            '--execute',
            '--to', 'notebook',
            '--output', tmp_file.name + '_executed',
            tmp_file.name
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if result.returncode == 0:
                return {'smoke_test': 'PASS', 'output': result.stdout}
            else:
                return {'smoke_test': 'FAIL', 'error': result.stderr}
        except subprocess.TimeoutExpired:
            return {'smoke_test': 'TIMEOUT', 'error': 'Execution timeout (60s)'}
        except Exception as e:
            return {'smoke_test': 'ERROR', 'error': str(e)}
```

### Phase 3: 総合評価とレポート生成

**評価スコア算出**:
```python
def calculate_quality_score(evaluation_results):
    """7軸評価の結果から総合品質スコアを算出（0-100）"""

    weights = {
        'syntax': 25,       # 構文エラーは致命的
        'dependencies': 20, # 依存関係問題も致命的
        'paths': 15,        # パス問題は修正可能
        'style': 5,         # スタイルは品質向上要因
        'reproducibility': 20,  # 再現性は重要
        'smoke_test': 10,   # 実行テストも重要
        'compatibility': 5  # 互換性は追加評価
    }

    scores = {}

    # 各軸のスコア計算
    scores['syntax'] = 0 if evaluation_results['syntax_errors'] else 25
    scores['dependencies'] = 20 * (1 - len(evaluation_results['missing_deps']) / max(evaluation_results['total_imports'], 1))
    scores['paths'] = 15 * (1 - len(evaluation_results['path_issues']) / max(len(evaluation_results['referenced_paths']), 1))
    scores['style'] = evaluation_results['style_score'] * 5 / 100
    scores['reproducibility'] = 20 if not evaluation_results['execution_issues'] else 10
    scores['smoke_test'] = 10 if evaluation_results['smoke_test'] == 'PASS' else 0
    scores['compatibility'] = 5  # 基本値（詳細評価は別途）

    total_score = sum(scores.values())
    return min(100, max(0, total_score))
```

## 評価レポートフォーマット

```markdown
# Code Evaluation Report

## 📊 総合評価

**総合スコア**: 85/100 ⭐⭐⭐⭐☆
**実行推奨度**: 🟢 RECOMMENDED

| 評価軸 | スコア | ステータス | 詳細 |
|-------|--------|-----------|------|
| 構文チェック | 25/25 | ✅ PASS | エラーなし |
| 依存関係 | 18/20 | ⚠️ MINOR | 2パッケージ要追加 |
| パス参照 | 12/15 | ⚠️ MINOR | 3ファイル要確認 |
| スタイル | 4/5 | ✅ GOOD | PEP8準拠率 92% |
| 再現性 | 20/20 | ✅ PASS | 実行順序問題なし |
| スモーク実行 | 10/10 | ✅ PASS | 基本動作確認済み |
| 環境互換性 | 4/5 | ✅ GOOD | 主要環境対応 |

## 🔧 修正が必要な問題

### Critical Issues (実行前に必修)
*問題なし*

### Major Issues (修正推奨)

**1. Missing Dependencies**
- `adjustText`: Volcano plotラベル調整用
- `gseapy`: パスウェイ解析用
```bash
pip install adjustText gseapy
```

**2. Path References**
- `data/raw/sample_info.csv`: サンプル情報ファイルが見つからない
- `results/` ディレクトリ: 自動作成されるが事前確認推奨

### Minor Issues (改善提案)

**1. Code Style**
- Line 45: Line too long (88 > 79 characters)
- Line 67: Missing whitespace after ','

## 🚀 実行準備チェックリスト

### 環境準備
- [ ] Python 3.8+ 環境
- [ ] `environment.yml` からconda環境作成
- [ ] 追加パッケージインストール: `adjustText`, `gseapy`

### データ準備
- [ ] `data/raw/` ディレクトリ作成
- [ ] mzMLファイルの配置確認
- [ ] sample_info.csv の作成

### 実行環境
- [ ] 8GB以上のRAM
- [ ] 約30分の実行時間確保
- [ ] ディスク空容量 5GB以上

## 🎯 実行成功予測

**予測実行成功率**: 92%

**想定される問題**:
1. 依存パッケージ未インストール (確率: 40%)
2. データファイル配置ミス (確率: 30%)
3. メモリ不足 (確率: 10%)

**推奨アクション**:
1. code-debugでの段階的実行確認
2. 依存関係の事前解決
3. データファイルの存在確認
```

## code-debug との連携

評価完了後、以下の情報をcode-debugに引き渡す：

```json
{
  "evaluation_status": "completed",
  "overall_score": 85,
  "execution_readiness": "READY_WITH_MINOR_FIXES",
  "critical_issues": [],
  "major_issues": [
    {
      "type": "missing_dependency",
      "packages": ["adjustText", "gseapy"],
      "fix_command": "pip install adjustText gseapy"
    }
  ],
  "minor_issues": [
    {"type": "style", "count": 2, "severity": "low"}
  ],
  "predicted_success_rate": 0.92,
  "recommended_debug_strategy": "incremental_execution",
  "checkpoint_cells": [5, 10, 15, 20],
  "high_risk_cells": []
}
```

## 高速評価モード

**軽量評価**（30秒以内）:
- 構文チェックのみ
- import文の存在確認のみ
- パス参照の基本チェックのみ

**完全評価**（5-10分）:
- 全7軸の詳細評価
- スモーク実行テスト
- 依存関係の詳細解析

評価完了後は**code-debug**で実際のJupyter Notebook実行・エラー修正・動作確認を行う。