# -*- coding: utf-8 -*-
import numpy as np
from scipy.interpolate import interp1d


class Force(object):
    def __init__(self, area, engine_data, T, density, Cd0=0.34):
        self.area = area
        self.engine_data = engine_data
        self.T = T
        self.density = density
        self.Cd0 = Cd0

        self.burn_time = np.max(engine_data[:, 0])
        self.edat = np.r_[self.engine_data, np.array([[self.T, 0.0]])]

    # 推力データを線形補間し計算時間分の1次元配列時系列データに整形
    def thrust(self):
        ft_t = interp1d(self.edat[:, 0], self.edat[:, 1])
        return ft_t

    def drag(self, alpha, air_velocity):
        Cd = self.Cd0 * (0.012 * (np.degrees(alpha))**2 + 1.0)
        V = np.linalg.norm(air_velocity)
        return 0.5 * Cd * self.density * self.area * V * air_velocity
