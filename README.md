## Overview
Xenium data is used. Sample data is below (https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer).
The shortest vector from one gene to another gene is extracted and projected onto polar coordinates. Then, PCoA is used 
to plot peripheral genes, internally expressed genes, and extensively expressed genes.

## install

```
git clone https://github.com/RyosukeNomural/Xenium.git
```

## Command Description


<br>

### 日本語版
## 概要
Xeniumデータを使用。サンプルデータは以下(https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer)。
ある遺伝子からみた他の遺伝子までの最短となるベクトルを抽出し、極座標上に投影する。その後PCoAによって、周辺遺伝子(peripheral)、内部発現遺伝子(internal)、広範囲発現遺伝子(ubiquitous)を図示する。

## インストール

```
git clone https://github.com/RyosukeNomural/Xenium.git
```

## コマンドの説明
- lib/preprocessing: csvファイルから本解析用にデータを加工。
- lib/mesh: 対象区間を10 μm ✕ 10 μm のメッシュに分割。解析対象遺伝子がどのメッシュに位置するかをマーキングし、各メッシュからどの解析対象遺伝子のメッシュまでが最短かを計算。
- lib/calculation: 遺伝子間の距離を計算。
- 
