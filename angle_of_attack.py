# -*- coding: utf-8 -*-
import numpy as np


# 迎角
def aoa(air_velocity):
    # 機軸座標系でz軸と対気速度ベクトルのなす角
    e = np.array([0.0, 0.0, 1.0])
    V = np.linalg.norm(air_velocity)
    if V == 0.0:
        return 0.0
    else:
        return np.arccos(np.dot(-air_velocity / V, e))
