i#! /usr/bin/env python

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

# DBSCANによって外れ値をフィルタリングする関数
def dbscan_filter(df, all_name_list: list, x_min, x_max, y_min, y_max, epithelial_gene, save_root):
    '''
    Parameters:
        df: 特定に区間に削ったdf
        all_name_list: dfに含まれる遺伝子名 dbscanを一つずつの遺伝子に適応する際に必要
        epithelial_gene : 始点に使う遺伝子名(大きく削る)
        save_root : 始点の塊を図示する場所

    Returns:
        None
        ただし、完成したcsvを保存する
    '''
    df_complete = pd.DataFrame()
    
    for name in all_name_list:
        if name == epithelial_gene:
            # 指定した遺伝子名のx,y座標を取り出す(dbscanでは数値のdfしか扱えないから)
            df_locations = df[df["feature_name"] == name].loc[:, 'x_location':'y_location']
            # 指定した遺伝子名の"transcript_id", "cell_id", "feature_name"を取り出す
            df_left = df[df["feature_name"] == name].loc[:, 'transcript_id':'feature_name']
            eps = 20
            N = 20
            dbscan = DBSCAN(eps=eps, min_samples=N)
            cluster_pred = dbscan.fit_predict(df_locations)
            df_locations['cluster'] = cluster_pred
            df_filtered = df_locations[(df_locations['cluster'] != -1)]

            # "transcript_id", "cell_id", "feature_name" と　'x_location', 'y_location',　"cluster"を連結する
            # clusterが形成されなかったものはデータフレーム自体がないので注意
            if not df_filtered.empty:
                df_ = pd.concat([df_left, df_filtered], axis=1)
                df_ = df_[pd.notnull(df_['x_location'])]
                df_complete = pd.concat([df_complete, df_])

        else:
            # 指定した遺伝子名のx,y座標を取り出す(dbscanでは数値のdfしか扱えないから)
            df_locations = df[df["feature_name"] == name].loc[:, 'x_location':'y_location']
            # 指定した遺伝子名の"transcript_id", "cell_id", "feature_name"を取り出す
            df_left = df[df["feature_name"] == name].loc[:, 'transcript_id':'feature_name']
            dbscan = DBSCAN(eps=300, min_samples=5)
            cluster_pred = dbscan.fit_predict(df_locations)
            df_locations['cluster'] = cluster_pred
            df_filtered = df_locations[(df_locations['cluster'] != -1)]

            # "transcript_id", "cell_id", "feature_name" と　'x_location', 'y_location',　"cluster"を連結する
            # clusterが形成されなかったものはデータフレーム自体がないので注意
            if not df_filtered.empty:
                df_ = pd.concat([df_left, df_filtered], axis=1)
                df_ = df_[pd.notnull(df_['x_location'])]
                df_complete = pd.concat([df_complete, df_])        

    # EPCAMの分布をplot        
    epcam_data = df_complete[df_complete['feature_name'] == epithelial_gene]

    # x_locationとy_locationの値を取得
    x_values = epcam_data['x_location']
    y_values = epcam_data['y_location']

    # 散布図をプロット
    
    plt.xlim(int(x_min), int(x_max))
    plt.ylim(int(y_min), int(y_max))  
    plt.axis('equal')
    plt.scatter(x_values, y_values, label=epithelial_gene, color='blue')


    # グラフにラベルやタイトルを追加
    plt.xlabel('x_location')
    plt.ylabel('y_location')
    plt.title(f'Scatter Plot of {epithelial_gene}')
    plt.savefig(f"{save_root}/{epithelial_gene}_scatter_plot_eps_{eps}_min_{N}")
    plt.close()
    print(df_complete.head(20))
                
    return df_complete

