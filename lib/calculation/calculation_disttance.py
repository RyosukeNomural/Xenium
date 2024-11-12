import numpy as np
# 距離を計算する関数
def calc_distance(x:float, y:float) -> float:
    return np.sqrt(x*x + y*y)
