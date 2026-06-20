---
name: code-implement
description: コード実装スキル - paper-selectで選定した論文のdry解析パイプラインを、仮想環境構築からコード実装まで、コピペで実行できる形で書く。無料かつ商用利用可能なツールのみ使用。有料ツールは無料代替で置き換える。「コード実装」「コードを書いて」「実装して」「解析コード」「パイプライン実装」「dry解析」「コーディング」等でトリガーする。
---

# Code Implementation Skill

論文のdry解析パイプラインを実装し、Jupyter Notebook形式で実行可能なコードを生成する。

## 実装ポリシー

### 1. 無料ツール優先ポリシー

**使用可能ツール（優先順位順）**:
```
1. Python標準ライブラリ（pandas, numpy, scipy, matplotlib, seaborn）
2. 無料・オープンソース専用ツール（sage-proteomics, pymzml, biopython）
3. アカデミック無料ツール（OpenMS, MaxQuant）
4. フリー版あり商用ツール（フリー機能のみ使用）
```

**置き換えルール**:
```
有料ツール → 無料代替案
- Progenesis QI → OpenMS Feature Detection
- Mascot → sage-proteomics / MSFragger
- PEAKS Studio → OpenMS + PyProphet
- Proteome Discoverer → sage-proteomics + MSstats
- Scaffold → MSstats + custom scripts
```

### 2. 実行可能性優先

**コード品質基準**:
- ✅ **コピペ実行可能**: 環境構築からデータ処理まで全て記述
- ✅ **依存関係明記**: pip/conda install コマンド完備
- ✅ **エラーハンドリング**: ファイル不在・ライブラリ未インストール対応
- ✅ **進行状況表示**: 長時間処理での進捗バー・ログ出力
- ✅ **中間ファイル保存**: 各ステップの結果をCSV/図として保存

## 実装ワークフロー

### Phase 1: 論文解析とツール選定

```
1. 論文のMethods sectionを詳細解析
   └→ 使用ツール・パラメータ・データフロー抽出
2. 有料ツールの無料代替案検索
   └→ 同等機能を持つオープンソースツールの特定
3. 実装可能性の評価
   └→ Python環境での実現可能性評価
4. データフロー設計
   └→ input → processing → output の明確化
```

### Phase 2: 環境設定コード生成

**environment.yml の自動生成**:
```yaml
# プロテオミクス解析環境のテンプレート
name: proteomics-analysis
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - numpy>=1.24
  - pandas>=2.0
  - scipy>=1.11
  - scikit-learn>=1.3
  - matplotlib>=3.7
  - seaborn>=0.13
  - biopython>=1.81
  - pymzml>=2.5
  - openms>=3.0  # 論文に応じて調整
  - pip
  - pip:
    - sage-proteomics>=0.14
    - gseapy>=1.0
    - adjustText>=0.8
```

### Phase 3: Jupyter Notebook 生成

**Notebook構造（標準テンプレート）**:

```python
# Cell 1: 環境確認とライブラリインポート
import sys, os, subprocess
import numpy as np
import pandas as pd
from pathlib import Path

# 必須ライブラリの存在確認
required_libs = ['numpy', 'pandas', 'scipy', 'matplotlib']
missing_libs = []
for lib in required_libs:
    try:
        __import__(lib)
        print(f"✅ {lib} available")
    except ImportError:
        missing_libs.append(lib)
        print(f"❌ {lib} missing")

if missing_libs:
    print(f"Missing libraries: {missing_libs}")
    print("Run: pip install " + " ".join(missing_libs))
    sys.exit(1)

# Cell 2: パス設定とディレクトリ作成
PROJECT_ROOT = Path(".").resolve()
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

# ディレクトリの自動作成
for dir_path in [DATA_DIR, RESULTS_DIR, FIGURES_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)
    print(f"📁 {dir_path} ready")

# Cell 3〜N: 論文の各Stepに対応するコード
# - データダウンロード
# - 前処理
# - 解析実行
# - 結果可視化
# - 統計解析
```

### Phase 4: 論文再現コード生成

**実装原則**:

1. **段階的実装**: 論文のMethodsを段階ごとに分割
2. **中間結果保存**: 各段階でCSV/pickle保存
3. **結果検証**: 論文のFigureと比較可能な図を生成
4. **パラメータ記録**: 使用したパラメータをJSON保存

**コード生成テンプレート**:

```python
# === Step X: [論文のステップ名] ===

def step_x_analysis(input_data, params):
    """
    論文のStep Xを再現する関数

    Parameters
    ----------
    input_data : DataFrame
        前段階からの入力データ
    params : dict
        解析パラメータ（論文記載値）

    Returns
    -------
    output_data : DataFrame
        解析結果データ
    figures : list
        生成された図のパス一覧
    """

    print(f"🔬 Step X: [ステップ名] 開始")
    print(f"📊 Input shape: {input_data.shape}")

    # --- 実際の解析処理 ---
    try:
        # 論文のアルゴリズムを実装
        output_data = process_data(input_data, params)

        # 中間結果の保存
        output_path = RESULTS_DIR / f"step_x_output.csv"
        output_data.to_csv(output_path, index=False)
        print(f"💾 Results saved: {output_path}")

        # 図の生成（論文のFigureを再現）
        fig_path = generate_figure(output_data)
        print(f"📈 Figure saved: {fig_path}")

        return output_data, [fig_path]

    except Exception as e:
        print(f"❌ Error in Step X: {e}")
        print("🔧 Troubleshooting suggestions:")
        print("   - Check input data format")
        print("   - Verify parameter values")
        raise

# 実行
step_params = {
    'parameter1': 'value_from_paper',
    'parameter2': 123.45
}

result_data, figures = step_x_analysis(input_data, step_params)
print(f"✅ Step X completed: {result_data.shape}")
```

## エラーハンドリング戦略

### 1. 実行前チェック

```python
# 依存関係チェック関数
def check_dependencies():
    """実行前に必要な環境をチェック"""

    checks = {
        'python_version': sys.version_info >= (3, 8),
        'required_libs': all(__import__(lib) for lib in REQUIRED_LIBS),
        'data_files': all(Path(f).exists() for f in REQUIRED_FILES),
        'disk_space': check_disk_space() > MIN_DISK_SPACE_GB
    }

    for check_name, passed in checks.items():
        status = "✅" if passed else "❌"
        print(f"{status} {check_name}")

    if not all(checks.values()):
        raise EnvironmentError("Environment check failed")
```

### 2. 段階的実行とロールバック

```python
# 実行状態の保存・復元
class ExecutionStateManager:
    """実行状態を管理し、エラー時のロールバックを可能にする"""

    def __init__(self, checkpoint_dir):
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(exist_ok=True)

    def save_checkpoint(self, step_name, data):
        """中間結果のチェックポイント保存"""
        checkpoint_path = self.checkpoint_dir / f"{step_name}.pkl"
        pd.to_pickle(data, checkpoint_path)
        print(f"💾 Checkpoint saved: {step_name}")

    def load_checkpoint(self, step_name):
        """チェックポイントから復元"""
        checkpoint_path = self.checkpoint_dir / f"{step_name}.pkl"
        if checkpoint_path.exists():
            data = pd.read_pickle(checkpoint_path)
            print(f"🔄 Restored from checkpoint: {step_name}")
            return data
        else:
            raise FileNotFoundError(f"No checkpoint found: {step_name}")
```

## code-evaluator との連携

実装完了後、code-evaluatorに以下の情報を引き渡す：

```json
{
  "implementation_status": "completed",
  "notebook_path": "/path/to/generated_notebook.ipynb",
  "environment_file": "/path/to/environment.yml",
  "required_data": ["file1.mzML", "file2.fasta"],
  "estimated_runtime": "15-30 minutes",
  "memory_requirement": "8GB RAM",
  "critical_dependencies": ["sage-proteomics>=0.14", "openms>=3.0"],
  "known_limitations": ["Requires GPU for >1000 samples", "Windows compatibility untested"]
}
```

## 出力仕様

**生成ファイル一覧**:
- `notebook_[paper_id].ipynb`: メインのJupyter Notebook
- `environment.yml`: conda環境定義
- `requirements.txt`: pip依存関係
- `README.md`: 実行手順書
- `config/`: パラメータ設定ファイル群
- `scripts/`: ユーティリティスクリプト群

**品質メトリクス**:
- ✅ コードカバレッジ: 論文のMethods sectionの90%以上をカバー
- ✅ 実行成功率: 依存関係満たした環境で100%実行成功
- ✅ 再現性: 論文のFigureと視覚的に一致する結果を生成
- ✅ ドキュメント: 各ステップに十分な説明コメント

実装完了後は**code-evaluator**で静的解析・動的評価を実施し、その後**code-debug**でJupyter Notebook上での実際のデバッグ・動作確認を行う。