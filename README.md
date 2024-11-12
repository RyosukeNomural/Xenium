## Overview
Xenium data is used. Sample data is below (https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer).
The shortest vector from one gene to another gene is extracted and projected onto polar coordinates. Then, PCoA is used 
to plot peripheral genes, internally expressed genes, and extensively expressed genes.

## install

```
git clone https://github.com/RyosukeNomural/Xenium.git
```

## Sample data
- lib/dataset/breast_sample/transcripts.csv: csv file before processing. Coordinates and expression information are stored.
- lib/dataset/breast_sample/transcripts_filter.csv: The csv file after processing. Use this file. Extract only necessary information.
Created with the following command.

```
cut --complement -f 3,5,7 -d ',' transcripts.csv > transcripts_filter.csv
```

  

## Command Description
- lib/preprocessing: Process data from csv files for this analysis. Low-expressed genes were filtered and clusters were extracted by DBSCAN.
- lib/mesh: Divides the interval into 10 μm✕ 10 μm meshes. It marks the mesh on which the genes to be analyzed are located and calculates the shortest distance from each mesh to the mesh of the gene to be analyzed.
- lib/calculation: Calculates the distance between genes.
- lib/similarity: Calculate similarity between polar coordinates.
- lib/pcoa: Project the table computed from the similarity onto a 2D space with dimensionality reduction.
- lib/result: where the results are stored. A part of polar plots of sample data and PCoA results (gene names are not shown) are stored.

<br>

### 日本語版
## 概要
Xeniumデータを使用。サンプルデータは以下(https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer)。
ある遺伝子からみた他の遺伝子までの最短となるベクトルを抽出し、極座標上に投影する。その後PCoAによって、周辺遺伝子(peripheral)、内部発現遺伝子(internal)、広範囲発現遺伝子(ubiquitous)を図示する。

## インストール

```
git clone https://github.com/RyosukeNomural/Xenium.git
```

## サンプルデータ
- lib/dataset/breast_sample/transcripts.csv: 加工前のcsvファイル。
- lib/dataset/breast_sample/transcripts_filter.csv: 加工後のcsvファイル。こちらを使用。必要な情報のみを抽出。
以下のコマンドで作成。

```
cut --complement -f 3,5,7 -d ',' transcripts.csv > transcripts_filter.csv
```

## コマンドの説明
- lib/preprocessing: csvファイルから本解析用にデータを加工。低発現遺伝子をフィルターし、DBSCANによりクラスターを抽出。
- lib/mesh: 対象区間を10 μm ✕ 10 μm のメッシュに分割。解析対象遺伝子がどのメッシュに位置するかをマーキングし、各メッシュからどの解析対象遺伝子のメッシュまでが最短かを計算。
- lib/calculation: 遺伝子間の距離を計算。
- lib/similarity: 極座標同士の類似度を計算。
- lib/pcoa: 類似度から計算したテーブルを次元削減し、2次元空間上に投影。
- lib/result: 結果の保存先。サンプルデータの極座標プロットの一部、およびPCoAの結果(遺伝子名は伏せてある)を保存してある。
