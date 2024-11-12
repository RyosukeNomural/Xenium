import numpy as np
from scipy.spatial import distance

def js_divergence(p: list, q:list)-> float:
    """
    p, q: 2つの確率分布（ヒストグラム）をnpのlistにしたもの
    """
    epsilon = 1e-10  # ゼロ除算を防ぐための微小な値
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)

    # ゼロ除算を避けるために微小な値を追加しておく
    p = p + epsilon
    q = q + epsilon

    # PとQの平均分布を計算
    m = 0.5 * (p + q)

    # KLダイバージェンスを計算
    kl_pm = np.sum(p * np.log(p / m))
    kl_qm = np.sum(q * np.log(q / m))

    # JSダイバージェンスを計算
    js = 0.5 * (kl_pm + kl_qm)
    return js

# braycurtis係数を全て計算する関数
def braycurtis(p: list, q: list) -> float:
    """
    p, q: 2つの確率分布（ヒストグラム）をnpのlistにしたもの
    """
    braycurtis_similarity_coefficient = distance.braycurtis(p,q)
    return braycurtis_similarity_coefficient
