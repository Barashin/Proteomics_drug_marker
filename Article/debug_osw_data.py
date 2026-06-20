#\!/usr/bin/env python3
import sqlite3, pandas as pd, os

def inspect_osw_schema(osw_file):
    conn = sqlite3.connect(osw_file)
    tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", conn)
    print(f"Tables in {os.path.basename(osw_file)}: {tables['name'].tolist()}")
    for table in ['FEATURE', 'TRANSITION', 'PRECURSOR']:
        try:
            schema = pd.read_sql_query(f"PRAGMA table_info({table})", conn)
            print(f"{table}: {[row['name'] for _, row in schema.iterrows()]}")
        except: pass
    conn.close()

inspect_osw_schema('results/openms_output/openswath/CRC01-N.osw')
