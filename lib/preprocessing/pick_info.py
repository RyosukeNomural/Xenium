#! /usr/bin/env python

from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import re
import os
from datetime import datetime

TRANSCRIPT_ID = 0
CELL_ID = 1
PROTEIN_NAME = 2
X_LOCATION = 3
Y_LOCATION = 4
Z_LOCATION = 5
QV = 6

## 特定の区間を切り出して、qvのフィルタリングをする関数
# txtファイルに書き出してパスを出力
def pick_location(FILE, x_min, x_max, y_min, y_max, panel):
    cut_file = ""
    # ファイル名に含む文字列
    area = f'{x_min}_{x_max}and{y_min}_{y_max}'
    print(area)
    x_min, x_max, y_min, y_max = float(x_min), float(x_max), float(y_min), float(y_max)

    with open(FILE, "r") as f:
        header = f.readline()
        for line in f:
            elems = line.strip('\n').split(',')
            protein_name = elems[PROTEIN_NAME].strip('"')

            if protein_name.startswith('NegControlCodeword') or protein_name.startswith('NegControlProbe') or protein_name.startswith('BLANK'):
                continue

            # qv20以上のみを抽出
            qv = float(elems[QV])
            if qv < 20:
                continue    

            # 範囲の限定
            x_location, y_location = float(elems[X_LOCATION]), float(elems[Y_LOCATION])
            if not (x_min <= x_location <= x_max and y_min <= y_location <= y_max):
                continue

            cut_file += line

        current_date = datetime.now()
        formatted_date = current_date.strftime("%Y%m%d%H%M%S")

        save_root = f"/home/nomu/work/Xenium/3_calc_min/lib/dataset/{panel}/{formatted_date}_{int(x_min)}_{int(x_max)}_{int(y_min)}_{int(y_max)}"
        if not os.path.exists(save_root): os.makedirs(save_root)
        
        file_name = f"{save_root}/{area}.txt"

    with open(file_name, "w") as f:
        f.write(cut_file)
    
    return file_name, save_root


## 遺伝子名順にソートする関数
def sort_file(file_name):
    df = pd.read_csv(file_name, names = ["transcript_id", "cell_id", "feature_name", "x_location", "y_location", "z_location", "qv"])
    # qvの値も消す
    del df['qv']
    
    # 'name', 'x_location'の順でデータフレームをソート
    df.sort_values(by=["feature_name", "x_location"], inplace=True)
    #print(df.head(10))
    return df

# データフレーム型から全ての遺伝子名を取得する関数
def get_gene_name(df)-> None:
    # df = pd.read_csv(FILE, names = ['name', 'x_location', 'y_location'])
    all_name_list = df.iloc[0:,2].unique().tolist()
    #print(all_name_list)
    print(f"all_name_list: {all_name_list}")
    return all_name_list
