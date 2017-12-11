# -*- coding: utf-8 -*-
import numpy as np
from scipy.interpolate import interp1d


class Mass(object):
    def __init__(self, m_af, I_af, CG_a, len_a, inertia_z0, inertia_zT, me_total, me_prop, len_e, d_e, burn_time, T):
        self.m_af = m_af
        self.I_af = I_af
        self.CG_a = CG_a
        self.len_a = len_a
        self.inertia_z0 = inertia_z0
        self.inertia_zT = inertia_zT
        self.me_total = me_total
        self.me_prop = me_prop
        self.len_e = len_e
        self.d_e = d_e
        self.burn_time = burn_time
        self.T = T

        self.CG_e = self.len_a - self.len_e * 0.5

        self.interp_time = np.array([0., self.burn_time, self.burn_time, self.T])
        self.M = np.array([self.m_af + self.me_total,
                           self.m_af + self.me_total - self.me_prop,
                           self.m_af + self.me_total - self.me_prop,
                           self.m_af + self.me_total - self.me_prop])
        self.Me = np.array([self.me_total,
                            self.me_total - self.me_prop,
                            self.me_total - self.me_prop,
                            self.me_total - self.me_prop])
        self.Me_dot = np.array([-self.me_total / self.burn_time,
                                -self.me_total / self.burn_time,
                                0.,
                                0., ])
        self.Inertia_z = np.array([self.inertia_z0,
                                   self.inertia_zT,
                                   self.inertia_zT,
                                   self.inertia_zT])
        self.Inertia_z_dot = np.array([(self.inertia_zT - self.inertia_z0) / self.burn_time,
                                       (self.inertia_zT - self.inertia_z0) / self.burn_time,
                                       0.,
                                       0.])

    # 全質量
    def mass(self):
        return interp1d(self.interp_time, self.M)

    # engine質量
    def me_t(self):
        return interp1d(self.interp_time, self.Me)

    # engine質量の時間変化量
    def me_dot(self):
        return interp1d(self.interp_time, self.Me_dot)

    # 全CG
    def CG(self):
        numerator = self.CG_a * self.m_af + self.CG_e * self.Me
        denominator = self.m_af + self.Me
        return interp1d(self.interp_time, numerator / denominator)

    # 全CGの時間変化量
    def CG_dot(self):
        CG_dot = self.Me_dot * (self.CG_e - self.CG_a) * self.m_af / ((self.m_af + self.Me)**2)
        return interp1d(self.interp_time, CG_dot)

    # engineのx, y慣性モーメント
    def iexg(self):
        ie = self.Me * (self.d_e * self.d_e / 16.0 + self.len_e * self.len_e / 12.0)
        return interp1d(self.interp_time, ie)

    # 全機体のx, y慣性モーメント
    def inertia(self):
        numerator = self.CG_a * self.m_af + self.CG_e * self.Me
        denominator = self.m_af + self.Me
        CG = numerator / denominator
        ie = self.Me * (self.d_e * self.d_e / 16.0 + self.len_e * self.len_e / 12.0)
        iner = self.I_af + self.m_af * ((self.CG_a - CG)**2) + ie + self.Me * ((self.CG_e - CG)**2)
        return interp1d(self.interp_time, iner)

    # z慣性モーメント
    def inertia_z(self):
        return interp1d(self.interp_time, self.Inertia_z)

    # x, y慣性モーメント時間微分
    def inertia_dot(self):
        numerator = self.CG_a * self.m_af + self.CG_a * self.Me
        denominator = self.m_af + self.Me
        CG = numerator / denominator
        CG_dot = self.Me_dot * (self.CG_e - self.CG_a) * self.m_af / ((self.m_af + self.Me) ** 2)
        ie_dot = self.Me_dot * (self.d_e * self.d_e / 16.0 + self.len_e * self.len_e / 12.0)
        iner_dot = self.m_af * 2 * (self.CG_a - CG) * (-CG_dot) \
                   + ie_dot + self.Me_dot * ((self.CG_e - CG)**2) + self.Me * 2 * (self.CG_e - CG) * (-CG_dot)
        return interp1d(self.interp_time, iner_dot)

    def inertia_z_dot(self):
        return interp1d(self.interp_time, self.Inertia_z_dot)
