import pandas as pd
import numpy as np

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
