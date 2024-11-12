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
import math
from lib.similarity.calc_divergence import js_divergence

def pickup_dist(gene_dict: dict, partner_gene: str) -> None:
    dists = []
    vectors = gene_dict[partner_gene]
    for vector in vectors:
        dists.append(math.sqrt(vector[0]**2 + vector[1]**2)
)
        
    #print(dists)
    return dists

# 始点から見たヒストグラムの値をリストに変換して保存する関数
def make_hist(gene_distances: list, start_gene: str, end_gene: str, dist_path: str, num_genes, selected_genes, colors) -> None:
    
    # 計算する範囲を指定
    range_min = 0
    range_max = 300

    # gene_distancesからinfを除外
    gene_distances = [x for x in gene_distances if x != float('inf')]
    # 範囲内のデータを取得
    gene_distances = [x for x in gene_distances if range_min <= x <= range_max]
    total_n = len(gene_distances)
    num_genes.append(total_n)

    if total_n > 50:
        selected_genes.append(end_gene) # 20以上のものだけを新たにリストへ
        colors.append(total_n)


        # CDH1から見た距離の分布をリストに
        from_hist, bin_edges, _ = plt.hist(gene_distances, bins=100, range=(0,300), density=True)
        #ここを変える~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if not os.path.exists(dist_path):
            os.makedirs(dist_path)
        save_path = f'{dist_path}/{end_gene}'

        plt.title(f'{end_gene}(n={total_n})')
        plt.xlabel('distance')
        plt.ylabel('frequency')
        plt.savefig(save_path)

        # リストに変換
        from_list = from_hist.tolist()
        plt.clf()
        plt.close()
        
        return from_list
    
# 距離の頻度が格納されたリストをある一定距離ごとに区切る
def split_hist(distribution: list):
    split_distances = []
    total_len = len(distribution) // 10
    for i in range(total_len):
        one_section = sum(distribution[i*10:i*10+10])
        split_distances.append(one_section)
    return(split_distances)

def dendrogram2newick(
    node: ClusterNode, parent_dist: float, leaf_names: List[str], newick: str = ""
) -> str:
    """Convert scipy dendrogram tree to newick format tree

    Args:
        node (ClusterNode): Tree node
        parent_dist (float): Parent distance
        leaf_names (List[str]): Leaf names
        newick (str): newick format string (Used in recursion)

    Returns:
        str: Newick format tree
    """
    if node.is_leaf():
        return f"{leaf_names[node.id]}:{(parent_dist - node.dist):.2f}{newick}"
    else:
        if len(newick) > 0:
            newick = f"):{(parent_dist - node.dist):.2f}{newick}"
        else:
            newick = ");"
        newick = dendrogram2newick(node.left, node.dist, leaf_names, newick)
        newick = dendrogram2newick(node.right, node.dist, leaf_names, f",{newick}")
        newick = f"({newick}"
        return newick

def pcoa_draw(root:str, standard_gene:str):
    """
    input: root(保存してあるファイルのroot)
           standard_gene(基準となる遺伝子)

    """

    # 遺伝子名とインデックスが一対一対応したディクショナリを作成
    root = root + "/"

    ## MAIN
    now = datetime.datetime.now()
    formatted_data = now.strftime("%Y%m%d")
    # save_path = f"/home/nomu/work/Xenium/3_calc_min/result/breast/20240612194030/2400_3300_2600_3300/fig_{formatted_data}_min_20/"
    save_path = root + f"fig_{formatted_data}_min_20/"
    if not os.path.exists(save_path): os.makedirs(save_path)

    f_path = root + "params.json"
    with open(f_path) as f:
        params = json.load(f)
        gene_list = params["gene_list"]
        print(f"全遺伝子数: {len(gene_list)}")

    # スパコン用
    distribution_path= root + "distribution/"

    # 全ての遺伝子の種類
    total_gene_len = len(gene_list)
    
    # cdh1からの全距離が格納されたファイルのパス 変える~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cdh1_file = root + standard_gene + "_dictionary.json"

    cdh1_open = open(cdh1_file, 'r')
    print("json file is opened.")
    cdh1_dict = json.load(cdh1_open)

    ## すべての遺伝子の距離の分布を格納
    all_gene_distribution = [] # すべての遺伝子の散布図をリスト化したものを格納する二重リスト
    split_all_gene_distribution = [] # 30ごとに区切ったすべての遺伝子の散布図をリスト化したものを格納する二重リスト
    sizes= []
    num_genes = []
    selected_genes = []
    for i in range(len(gene_list)):
        # 距離を格納するリスト
        partner_distances = []
        partner_gene = gene_list[i]

        # 1つの遺伝子について、cdh1からの全距離を抽出
        min_distances = pickup_dist(cdh1_dict, partner_gene)
        gene_distribution = make_hist(min_distances, standard_gene, partner_gene, distribution_path, num_genes, selected_genes, sizes)
        if gene_distribution:
            all_gene_distribution.append(gene_distribution)
            split_distance = split_hist(gene_distribution) # 30ごとに区切る
            split_all_gene_distribution.append(split_distance)
    print(f"all_genae_distribution: {len(all_gene_distribution)}")
    print("loading is finished!!!")

    #cdh1からの全距離はもう不要なので削除
    del cdh1_dict
    gc.collect()

    selected_genes = [item for item in selected_genes if 'CDH1' not in item]
    total_gene_len = len(selected_genes)
    print(f"total_gene_len: {total_gene_len}")

    ## ==================すべてのjsを計算=============================
    js_table = [[[0] for _ in range(total_gene_len)] for _ in range(total_gene_len)] # n×nの係数を入れる二次元行列
    for i in range(total_gene_len):
        for j in range(total_gene_len):
            standard_list = np.array(all_gene_distribution[i])
            partner_list = np.array(all_gene_distribution[j])
            js_ = js_divergence(standard_list, partner_list)
            js_table[i][j] = js_


    js_matrix = pd.DataFrame(js_table,columns=selected_genes, index=selected_genes)
    js_matrix = js_matrix.fillna(99999)
    print(js_matrix)
    js_matrix.to_csv(f"{save_path}/js_divergence_over_50_non_CDH1.csv")
    df = js_matrix

    upper_triangular_matrix = df.where(np.triu(np.ones(df.shape), k=1).astype(bool))
    upper_triangular_matrix = upper_triangular_matrix.fillna(0)
    # 行列形式に変換
    upper_triangular_matrix = upper_triangular_matrix.values

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
    plt.savefig(f"{save_path}/PCoA_Eigenvalue_bray_non_CDH1_over50")
    plt.show()
    plt.close()
