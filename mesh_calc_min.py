#!/usr/bin/env python

from lib.preprocessing.csv_to_dataframe import csv_to_dataflame
from lib.calculation.calculation_distance import calc_distance
from lib.make_file import make_file
from lib.save.save_to_json import save_params, save_distances
import sys
from collections import defaultdict
import numpy as np
import time
import gc
import json
import sys
import pandas as pd
import os
from datetime import datetime

TRANS_ID = 0 # 各mRNAのID
CELL_ID  = 1 # 各細胞のID
PROTEIN_NAME = 2 # mRNA名
X_LOCATION = 3 # x座標
Y_LOCATION = 4 # y座標


def mesh_calc_min_vector(x_min, x_max, y_min, y_max, epithelial_gene, panel): # sysで入力 
    gene_name_list, FILE, panel =  make_file(x_min, x_max, y_min, y_max, epithelial_gene, panel)  

    current_date = datetime.now()
    formatted_date = current_date.strftime("%Y%m%d%H%M%S")
    path_root = f"/home/nomu/work/Xenium/3_calc_min/result/min_vector/{panel}/{formatted_date}/"
    # パスが存在しない場合、新しく作成する
    if not os.path.exists(path_root):
        os.makedirs(path_root)
        print(f'Made {path_root} directory!')
    else:
        print(f'You already have {path_root} directory!')

    #### 4点の座標から正方形のメッシュを作成 ####
    MESH_SIZE = 10
    params=sys.argv
    x_left = x_min
    x_right = x_max
    y_bottom = y_min
    y_top = y_max

    x_len = x_right - x_left # length in x-direction
    y_len = y_top - y_bottom # length in y-direction
    x_num = int((x_len + MESH_SIZE) // MESH_SIZE)  # division number in x-direction(切り上げ)
    y_num = int((y_len + MESH_SIZE) // MESH_SIZE) # division number in y-direction

    mesh_numbers = []
    x_pointers, y_pointers = [], []
    mesh_lefts, mesh_rights, mesh_bottoms, mesh_tops = [], [], [], [] # メッシュのxの小さい方,大きい方,y座標の小さい方,大きいほうをしまう
    mesh_number = 0 # 何番目のメッシュに該当するか((0,0)→(0,1)→(0,2)の順番で1,2,3と振っていく)

    for i in range(y_num - 1):
        bottom = y_bottom + i * MESH_SIZE #  # y座標の小さい方を10ずつ検索
        top = bottom + MESH_SIZE # y座標の大きい方を10ずつ検索
        y_pointer = int(bottom // MESH_SIZE) # y座標の小さい方を10で割った値をポインターにする

        for j in range(x_num - 1):
            left = x_left + j * MESH_SIZE
            right = left + MESH_SIZE
            x_pointer = int(left // MESH_SIZE) # x座標の小さい方を10で割った値をポインターにする
            mesh_number +=1

            mesh_bottoms.append(bottom)
            mesh_tops.append(top)
            mesh_lefts.append(left)
            mesh_rights.append(right)
            mesh_numbers.append(mesh_number)
            x_pointers.append(x_pointer)
            y_pointers.append(y_pointer)
            #print(f"{mesh_number}:{left}_{right}_{bottom}_{top}")

    mesh_df = pd.DataFrame(
            data={'mesh_pointer': mesh_numbers, 'left': mesh_lefts, 'right': mesh_rights, 
                'bottom':mesh_bottoms, 'top': mesh_tops, 'x_pointer': x_pointers, 'y_pointer': y_pointers,
                'x_y_pointer': list(zip(x_pointers, y_pointers)) },
            columns=['mesh_pointer', 'left', 'right', 'bottom', 'top', 'x_pointer', 'y_pointer', 'x_y_pointer'],
            index=None
        )

    mesh_df.to_csv(f'{path_root}/mesh.csv', index=False)

    ### transcipts_filteredの読み込み
    gene_df = pd.read_csv(FILE)
    gene_df = gene_df.sort_values('x_location')
    gene_df.to_csv(f'{path_root}/gene_mesh.csv', index=False)

    ## 各x,y座標からどのポインターに対応するかを決定
    gene_df['x_y_pointer'] = list(zip(gene_df['x_location'].astype(int) // 10, gene_df['y_location'].astype(int) // 10))

    ### CDH1のみのdfを作成 
    cdh1_df = gene_df[gene_df["feature_name"] == "CDH1"]

    ### cdh1のポインターを抽出
    cdh1_pointers = cdh1_df["x_y_pointer"].unique()

    ### x_y_pointerごとにグループ化し、各グループをリストとして保存
    cdh1_grouped = cdh1_df.groupby('x_y_pointer')

    ### グループごとのDataFrameを辞書として保存
    cdh1_grouped_dfs_dict = {group: df for group, df in cdh1_grouped}

    ### gene_dfから全てのポインターを抽出(gene_pointers)
    gene_pointers = gene_df["x_y_pointer"].unique()

    ### gene_dfのポインターからcdh1のポインターまでの距離を計算し、{どこから: どこへ} が最小なのかのディクショナリ(from_gene_to_cdh1)を作成
    gene_names = sorted(gene_df["feature_name"].unique())
    print(f"gene_names")

    ### あるメッシュからどのCDH1のメッシュまでが最短かを計算
    from_gene_to_cdh1 = {}

    # gene_dfのポインターを抽出
    for gene_pointer in gene_pointers:
        min_distance = float('inf')
        keep_cdh1_pointer = ""

        # cdh1のポインターを抽出
        for cdh1_pointer in cdh1_grouped_dfs_dict.keys():
            # ポインターごとの距離を計算
            distance = np.square(gene_pointer[0] - cdh1_pointer[0]) + np.square(gene_pointer[1] - cdh1_pointer[1])
            # 最小値を取る時のポインターを保存
            if distance <= min_distance:
                keep_cdh1_pointer = cdh1_pointer
                min_distance = distance
        from_gene_to_cdh1[gene_pointer] = keep_cdh1_pointer


    ##### MAIN #####
    distance_dict = {} # 最小距離をすべて格納
    time_sta = time.time()

    ### gene_dfを遺伝子名ごとに参照し、新たなdfにする
    for start_gene_name in gene_names:
        mins = []
        if start_gene_name != "CDH1":
            start_gene_df = gene_df[gene_df['feature_name'] == start_gene_name]

            ### ポインターを読んで、どのCDH1のポインターを読めばいいかディクショナリから探す
            for i in range(start_gene_df.shape[0]):
                vector = []
                min_distance = float('inf')
                std_x, std_y = start_gene_df.iloc[i, X_LOCATION], start_gene_df.iloc[i, Y_LOCATION]

                start_gene_pointer = start_gene_df.iloc[i, 6]
                cdh1_target_pointer = from_gene_to_cdh1[start_gene_pointer]

                ## 必要箇所のCDH1
                filtered_epithelial_df = cdh1_grouped_dfs_dict[cdh1_target_pointer]

                # 最短距離計算
                for j in range(filtered_epithelial_df.shape[0]):
                    x_location, y_location = filtered_epithelial_df.iloc[j, X_LOCATION], filtered_epithelial_df.iloc[j, Y_LOCATION]
                    # 始点から終点までの距離
                    protein_distance = calc_distance(std_x - x_location, std_y - y_location)
                    # 最小距離の更新
                    if protein_distance <= min_distance:
                        min_distance = protein_distance
                        vector = [std_x - x_location, std_y - y_location]
                if vector != []:
                    mins.append(vector)
            print(start_gene_name)
            distance_dict[start_gene_name] = mins
            print(time.time() - time_sta)
    print("calculation is completed!")

    area = f"{x_min}_{x_max}_{y_min}_{y_max}"

    path_name = f'{path_root}/{area}/'
    # パスが存在しない場合、新しく作成する
    if not os.path.exists(path_name):
        os.makedirs(path_name)
        print('Made this new directory!')
    else:
        print('You already have this directory!')


    file_name = path_name + f'{epithelial_gene}_dictionary.json'
    # 距離をdictionaryで格納
    save_distances(file_name, distance_dict)
    # 遺伝子のリスト、保存先のroot、保存名をjsonファイルに格納
    save_params(gene_name_list, path_name, file_name)

    return gene_name_list, path_name, file_name


if __name__ == "__main__":
    #mesh_calc_min_vector(5720, 5920, 3140, 3410, "CDH1", "breast/case1_after")
    #mesh_calc_min_vector(2530, 2800, 270, 480, "CDH1", "lung")
    #mesh_calc_min_vector(4000, 5200, 200, 650, "KRT7", "/breast/case7_after/")
    mesh_calc_min_vector(3300, 3700, 3080, 3500, "CDH1", "breast/case1_after")


