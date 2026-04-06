#!/usr/bin/env python3
"""
Step 1: データ取得
================
Toyota et al. 2025 の公開データをダウンロードするスクリプトです。

【ダウンロードするもの】
  1. 論文の補足データ（Supplementary Tables S1〜S18）
     - Table S2 がメインの入力データ（タンパク質定量マトリクス）
  2. UniProt ヒトプロテオーム FASTA
     - OpenMS等のDIA解析ソフトでタンパク質同定に使用する配列データベース

【データの出典】
  - 論文: Toyota et al., Proteomes 2025; 13(3):38
  - DOI: https://doi.org/10.3390/proteomes13030038
  - ProteomeXchange ID: PXD058672
  - jPOST ID: JPST003422

入力: なし（インターネット接続が必要）
出力: data/raw/ 以下にファイルを保存
"""

# ============================================================
# ライブラリの読み込み
# ============================================================
import os        # ファイルパス操作・ディレクトリ作成に使用
import requests  # HTTPリクエスト（ファイルのダウンロード）に使用
import sys       # Pythonのバージョン確認等に使用

# ============================================================
# 設定（パスやURLの定義）
# ============================================================

# データ保存先ディレクトリ
# __file__ はこのスクリプト自身のパスを指す
# os.path.join() でOS（Mac/Linux/Windows）に依存しないパスを作成
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "raw")

# ディレクトリが存在しなければ作成（exist_ok=True で既に存在してもエラーにならない）
os.makedirs(DATA_DIR, exist_ok=True)

# 論文の補足データURL（MDPIの論文ページ）
BASE_URL = "https://www.mdpi.com/2227-7382/13/3/38"
# MDPI のリソースサーバーから直接ダウンロード（論文ページ経由だと403エラーになる場合がある）
SUPP_URL = "https://mdpi-res.com/d_attachment/proteomes/proteomes-13-00038/article_deploy/proteomes-13-00038-s001.zip"

# UniProt ヒトプロテオーム FASTA のダウンロードURL
# proteome:UP000005640 = ヒト（Homo sapiens）の参照プロテオーム
# compressed=true: gzip圧縮でダウンロード（サイズ削減）
FASTA_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    "compressed=true&format=fasta&query=proteome:UP000005640"
)


# ============================================================
# ダウンロード関数
# ============================================================
def download_file(url, filepath, description=""):
    """
    指定URLからファイルをダウンロードして保存する関数。

    引数:
        url (str): ダウンロード元のURL
        filepath (str): 保存先のファイルパス
        description (str): 表示用の説明文（省略可）

    戻り値:
        bool: ダウンロード成功なら True、失敗なら False
    """

    # 既にファイルが存在する場合はダウンロードをスキップ
    # （何度実行しても安全なように、冪等性を確保）
    if os.path.exists(filepath):
        print(f"  スキップ（既存）: {os.path.basename(filepath)}")
        return True

    print(f"  ダウンロード中: {description or os.path.basename(filepath)}")

    try:
        # requests.get() でHTTP GETリクエストを送信
        # timeout=60: 60秒以内に応答がなければタイムアウト
        # allow_redirects=True: リダイレクト（転送）を自動追従
        resp = requests.get(url, timeout=60, allow_redirects=True)

        # ステータスコードが200番台以外（エラー）なら例外を発生
        resp.raise_for_status()

        # ダウンロードしたデータをバイナリモード("wb")でファイルに書き込み
        with open(filepath, "wb") as f:
            f.write(resp.content)

        # ファイルサイズを MB 単位で表示
        size_mb = len(resp.content) / (1024 * 1024)
        print(f"  完了 ({size_mb:.1f} MB)")
        return True

    except requests.exceptions.Timeout:
        print(f"  エラー: タイムアウト（60秒以内に応答がありませんでした）")
        return False
    except requests.exceptions.HTTPError as e:
        print(f"  エラー: HTTPエラー ({e})")
        return False
    except Exception as e:
        print(f"  エラー: {e}")
        return False


# ============================================================
# メイン処理
# ============================================================
def main():
    print("=" * 60)
    print("Step 1: Toyota et al. 2025 データ取得")
    print("=" * 60)

    # ----------------------------------------------------------
    # [1/2] 論文の補足データ（Supplementary Tables）のダウンロード
    # ----------------------------------------------------------
    # Table S2 に DIA解析の出力と同等のタンパク質定量データが含まれている。
    # DIA解析ステップ（Step 03）をスキップしたい場合は、この Table S2 から解析を始められる。
    print("\n[1/2] 論文の補足データをダウンロード")
    print(f"  URL: {SUPP_URL}")

    supp_path = os.path.join(DATA_DIR, "proteomes-13-00038-s001.zip")
    success = download_file(SUPP_URL, supp_path, "Supplementary Tables (zip)")

    # ダウンロード成功＆ZIPファイルの場合は自動展開
    if success and supp_path.endswith(".zip"):
        import zipfile  # ZIP展開ライブラリ（必要な時だけインポート）
        try:
            # ZIPファイルを開いて、中身を data/raw/ に展開
            with zipfile.ZipFile(supp_path, "r") as z:
                z.extractall(DATA_DIR)
            print("  展開完了")
        except zipfile.BadZipFile:
            # MDPIサーバーからHTMLが返ってきた場合など
            print("  注意: ZIPファイルではない可能性があります。")
            print("  手動で論文ページからダウンロードしてください。")

    # ----------------------------------------------------------
    # [2/2] UniProt ヒトプロテオーム FASTA のダウンロード
    # ----------------------------------------------------------
    # FASTA ファイル: タンパク質の名前とアミノ酸配列が記録されたテキストファイル。
    # OpenMS等のDIA解析ソフトがペプチドを同定する際の「辞書」として使用する。
    # UniProt (https://www.uniprot.org) は世界最大のタンパク質配列データベース。
    print("\n[2/2] UniProt ヒトプロテオームFASTAをダウンロード")
    fasta_path = os.path.join(DATA_DIR, "human_proteome.fasta.gz")
    download_file(FASTA_URL, fasta_path, "UniProt UP000005640 (human proteome)")

    # ----------------------------------------------------------
    # 完了メッセージとファイル一覧の表示
    # ----------------------------------------------------------
    print("\n" + "=" * 60)
    print("データ取得完了！")
    print(f"保存先: {os.path.abspath(DATA_DIR)}")
    print()

    # ダウンロードしたファイルの一覧とサイズを表示
    print("ファイル一覧:")
    for f in sorted(os.listdir(DATA_DIR)):
        size = os.path.getsize(os.path.join(DATA_DIR, f))
        # サイズを読みやすい形式で表示（例: 1,234,567 bytes）
        print(f"  {f} ({size:,} bytes)")

    print()
    print("注意: 補足データが正しくダウンロードできない場合は、")
    print(f"  {BASE_URL} から手動でダウンロードしてください。")
    print("=" * 60)


# ============================================================
# スクリプトの実行エントリーポイント
# ============================================================
# このファイルを直接実行した場合のみ main() を呼び出す。
# 他のスクリプトから import した場合は実行されない。
if __name__ == "__main__":
    main()
