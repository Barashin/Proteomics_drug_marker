#!/usr/bin/env python3
"""
OSWファイルのスキーマ調査
"""

import sqlite3
import sys

def inspect_osw_schema(osw_file):
    """OSWファイルのスキーマを調査"""
    print(f"=== {osw_file} スキーマ調査 ===")

    try:
        conn = sqlite3.connect(osw_file)
        cursor = conn.cursor()

        # テーブル一覧
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [row[0] for row in cursor.fetchall()]
        print(f"テーブル数: {len(tables)}")
        print(f"テーブル: {', '.join(tables)}")

        # 主要テーブルの構造確認
        key_tables = ['FEATURE', 'TRANSITION', 'PRECURSOR', 'PEPTIDE', 'PROTEIN', 'RUN']

        for table in key_tables:
            if table in tables:
                print(f"\n--- {table} テーブル ---")
                cursor.execute(f"PRAGMA table_info({table});")
                columns = cursor.fetchall()
                for col in columns:
                    print(f"  {col[1]} ({col[2]})")

                # サンプルデータ確認
                cursor.execute(f"SELECT COUNT(*) FROM {table};")
                count = cursor.fetchone()[0]
                print(f"  レコード数: {count}")

                if count > 0:
                    cursor.execute(f"SELECT * FROM {table} LIMIT 3;")
                    sample_data = cursor.fetchall()
                    print(f"  サンプル: {len(sample_data)} rows")

        conn.close()

    except Exception as e:
        print(f"エラー: {e}")

if __name__ == "__main__":
    # 最初のOSWファイルを調査
    osw_file = "results/openms_output/openswath/CRC01-N.osw"
    inspect_osw_schema(osw_file)