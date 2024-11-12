#!/usr/bin/env python

import datetime
import json
import matplotlib.pyplot as plt
import sys
import numpy as np
import os
from scipy.spatial import distance
from scipy.signal import find_peaks
import seaborn as sns
import gc
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as shc
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import zscore
from skbio.stats.ordination import pcoa
import skbio
from mpl_toolkits.mplot3d import Axes3D
from sklearn.manifold import TSNE
from ete3 import Tree
import scipy.cluster.hierarchy as hc
from scipy.cluster.hierarchy import ClusterNode
from typing import List
from scipy.cluster.hierarchy import to_tree

def make_hist(vectors: list, start_gene: str, end_gene: str, dist_path: str, num_genes, selected_genes, colors) -> None:
    
    # 計算する範囲を指定
    range_min = 0
    range_max = 600

    # vectorsからinfを除外
    vectors = [v for v in vectors if all(np.isfinite(v))]
    
    # 範囲内のデータを取得
    vectors = [v for v in vectors if range_min <= np.sqrt(v[0]**2 + v[1]**2) <= range_max]
    total_n = len(vectors)
    num_genes.append(total_n)

    if total_n > 1:
        selected_genes.append(end_gene) # 50以上のものだけを新たにリストへ
        colors.append(total_n)

        ax = plt.subplot(111, projection="polar")

        # x, y から r, theta への変換
        angles = [np.arctan2(y, x) for x, y in vectors]
        radii = [np.sqrt(x**2 + y**2) for x, y in vectors]

        # ヒートマップ用のデータ作成
        heatmap, r_edges, theta_edges = np.histogram2d(radii, angles, bins=[20, 12], range=[[0, 600], [-np.pi, np.pi]])
        #[半径][角度]で入っている
        #print(r_edges) 
        #[  0.   5.  10.  15.  20.  25.  30.  35.  40.  45.  50.  55.  60.  65..... 300.]
        
        # ビンの個数を全体のデータ数で割り、割合を計算
        heatmap_percentage = heatmap / total_n
        #print(heatmap_percentage)
        #print("================================")

        
        # r, theta の中心を計算
        r_centers = (r_edges[:-1] + r_edges[1:]) / 2
        theta_centers = (theta_edges[:-1] + theta_edges[1:]) / 2

        # メッシュグリッドの作成
        R, Theta = np.meshgrid(r_centers, theta_centers)

        # ヒートマップの描画
        c = ax.pcolormesh(Theta, R, heatmap_percentage.T, shading='auto', cmap='viridis')

        # カラーバーを追加
        fig = plt.gcf()
        fig.colorbar(c)

        # CDH1から見た距離の分布をリストに   
        # 保存ディレクトリの作成
        if not os.path.exists(dist_path):
            os.makedirs(dist_path)
        save_path = os.path.join(dist_path, f'{end_gene}_dense_.png')

        # プロットを保存
        plt.title(f'{end_gene}(n={total_n})')
        plt.savefig(save_path, dpi=300)
        

        # プロットをクリアして閉じる
        plt.clf()
        plt.close()
        
        return heatmap_percentage

# 最小距離を抽出する関数
def pickup_vectors(gene_dict: dict, partner_gene: str) -> None:
    vectors = [] # ある遺伝子までの距離をすべて格納するリスト
    vectors = gene_dict[partner_gene]
    return vectors

def sum_matrix(matrix):
    return sum(sum(row) for row in matrix)

def min_matrix(list1, list2):
    if len(list1) != len(list2) or any(len(row1) != len(row2) for row1, row2 in zip(list1, list2)):
        raise ValueError("Input lists must have the same dimensions")

    result = [[min(list1[i][j], list2[i][j]) for j in range(len(list1[i]))] for i in range(len(list1))]
    total_sum = sum_matrix(result)

    return total_sum 

# NumPy配列から2次元リストへの変換
def numpy_array_to_list(arr):
    return arr.tolist() if isinstance(arr, np.ndarray) else arr

# PCoAを実行するコード
def calc_pcoa(root:str, standard_gene:str):
    ## MAIN
    now = datetime.datetime.now()
    formatted_data = now.strftime("%Y%m%d")
    save_path = f"{root}/fig_{formatted_data}/"
    
    
    f_path = root + "params.json"
    with open(f_path) as f:
        params = json.load(f)
        gene_list = params["gene_list"]
        gene_list.remove("CDH1")
        print(f"全遺伝子数: {len(gene_list)}")
        
    # 遺伝子名とインデックスが一対一対応したディクショナリを作成
    # cdh1からの全距離が格納されたファイルのパス 変える~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cdh1_file = root + standard_gene + "_dictionary.json"
    
    cdh1_open = open(cdh1_file, 'r')
    print("json file is opened.")
    cdh1_dict = json.load(cdh1_open)
    
    # スパコン用
    distribution_path= root + "distribution/"
    
    # 全ての遺伝子の種類
    total_gene_len = len(gene_list)
    
    ## すべての遺伝子の距離の分布を格納
    all_gene_distribution = [] # すべての遺伝子の散布図をリスト化したものを格納する二重リスト
    sizes= []
    num_genes = []
    selected_genes = []
    for i in range(len(gene_list)):
        # 距離を格納するリスト
        partner_distances = []
        partner_gene = gene_list[i]
    
        # 1つの遺伝子について、cdh1からの全距離を抽出
        min_vectors = pickup_vectors(cdh1_dict, partner_gene)
        heatmaps = make_hist(min_vectors, standard_gene, partner_gene, distribution_path, num_genes, selected_genes, sizes)
        if heatmaps is not None:
            all_gene_distribution.append(heatmaps)

    ## ==================すべてのjsを計算=============================
    js_table = [[[0] for _ in range(total_gene_len)] for _ in range(total_gene_len)] # n×nの係数を入れる二次元行列
    for i in range(total_gene_len):
        for j in range(total_gene_len):
            standard_list = numpy_array_to_list((all_gene_distribution[i]))
            partner_list = numpy_array_to_list((all_gene_distribution[j]))
            js_ =  min_matrix(standard_list, partner_list)
            js_table[i][j] = js_

    js_matrix = pd.DataFrame(js_table,columns=selected_genes, index=selected_genes)
    js_matrix = js_matrix.fillna(0)
    #print(js_matrix)
    df = js_matrix
    df = 1 - df # 対角化成分が0になるように変換
    np.fill_diagonal(df.values, 0)
    #print(df)
    
    upper_triangular_matrix = df.where(np.triu(np.ones(df.shape), k=1).astype(bool))
    upper_triangular_matrix = upper_triangular_matrix.fillna(0)
    # 行列形式に変換
    upper_triangular_matrix = upper_triangular_matrix.values
    print(upper_triangular_matrix)

    ## PCoA
    # 上三角行列からDistanceMatrixを作成
    triangle_values = df.values
    distance_matrix = skbio.DistanceMatrix(triangle_values)
    
    # PCoAを適用
    pcoa_results = pcoa(distance_matrix)
    print(pcoa_results)
    # PCoAの結果から座標を取得
    coordinates = pcoa_results.samples

    # 寄与率の取得
    explained_variance_ratio = pcoa_results.proportion_explained
    eigenvalues = explained_variance_ratio[:20]
    # 結果を表示
    print("Explained Variance Ratio:")
    print(explained_variance_ratio)
    
    # 座標軸の番号
    axis_numbers = range(1, 21)
    x_axis_numbers = ["PCoA" + str(i) for i in axis_numbers]
    
    # 棒グラフを描画
    plt.bar(x_axis_numbers, eigenvalues, color='skyblue')
    plt.ylabel('Eigenvalue')
    plt.xticks(rotation='vertical',fontsize=8)
    plt.subplots_adjust(left=0.13, right=0.99, bottom=0.18, top=0.99)
    plt.savefig(f"{root}/pcoa_eigenvalue")
    plt.show()
    plt.close()

    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    # 二次元プロット
    plt.scatter(coordinates.iloc[:, 0], coordinates.iloc[:, 1], s=80)
    fig_names = [""]
    # 各点にラベルを付ける
    for i, sample_id in enumerate(df.index):
        if sample_id in fig_names:
            plt.annotate(sample_id, (coordinates.iloc[i, 0], coordinates.iloc[i, 1]))

    # プロットの装飾
    plt.xlabel('PCoA1', fontsize=22)
    plt.xticks(fontsize=17)
    plt.ylabel('PCoA2', fontsize=22)
    plt.yticks(fontsize=17)
    plt.title('PCoA Plot')
    plt.show()
    plt.close()

    return coordinates

