#! /usr/bin/env python

from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import re
from sklearn.cluster import DBSCAN
import os
from datetime import datetime

# 3000以上ある遺伝子は3000をランダムに抽出して計算
def remove_house_keeping(df):
    '''
    Parameters:
        df: dbscanを行って外れ値を除去したdf
        file_name: txtファイルの名前 csvファイルのもとにする

    Returns:
        None
    '''
    # 同一遺伝子の上限数(5000)
    threshold = 5000

    # 遺伝子名ごとに行数をカウント
    value_counts = df['feature_name'].value_counts()
    #print(value_counts)

    house_keeping_name = []

    # 1列目の値が3000個以上の場合、ランダムに3000個のみを抽出
    # 3000個未満の場合はそのまま全て残す
    result_df = pd.DataFrame()
    for value, count in value_counts.items():
        if count >= threshold:
            sampled_rows = df[df['feature_name'] == value].sample(n=threshold, random_state=42)
            result_df = pd.concat([result_df, sampled_rows])
            house_keeping_name.append(value)
        else:
            result_df = pd.concat([result_df, df[df['feature_name'] == value]])

    house_keeping_name.sort()
    print(f'{house_keeping_name} are reduced to 3000!')
    return result_df

def remove_rare(df):
    # 同一遺伝子の下限数(5)
    threshold = 5
    # 遺伝子名ごとに行数をカウント
    value_counts = df['feature_name'].value_counts()

    for value, count in value_counts.items():
        print(f"name:{value} count:{count}")
        if count <= threshold:
            df = df[df['feature_name'] != value]

    # dfから条件に合致した行が削除された結果を返す
    return df

