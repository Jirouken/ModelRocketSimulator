# -*- coding: utf-8 -*-
import numpy as np
"""
クォータニオンをベクトルq = [q0, q1, q2, q3]で定義する．
共役クォータニオンq* = [q0, -q1, -q2, -q3]
クォータニオン積を行列表記すると，
qp = [[q0, -q1, -q2, -q3], [p0,
      [q1,  q0, -q3,  q2],  p1,
      [q2,  q3,  q0, -q1],  p2,
      [q3, -q2,  q1,  q0]]  p3]
ある座標系のベクトルr_0をq_01によって座標系を回転させると，
r_1 = q_01* r_1 q_01
"""


class Quaternion(object):
    # 3次元ベクトルuをクォータニオンに変換.
    def quaternion(self, u):
        return np.r_[0.0, u]

    # クォータニオンqの共役クォータニオン
    def coquat(self, q):
        return np.array([q[0], -q[1], -q[2], -q[3]])

    # クォータニオンの掛け算
    def product(self, q1, q2):
        q1 = np.array([[q1[0], -q1[1], -q1[2], -q1[3]],
                       [q1[1],  q1[0], -q1[3],  q1[2]],
                       [q1[2],  q1[3],  q1[0], -q1[1]],
                       [q1[3], -q1[2],  q1[1],  q1[0]]])
        return np.dot(q1, q2)

    # ベクトルrをクォータニオンqによって座標系を回転させたベクトルを返す.
    def rotation(self, r, q):
        r = self.quaternion(r)
        coq = self.coquat(q)
        r_ = self.product(coq, self.product(r, q))
        return np.array([r_[1], r_[2], r_[3]])

    # クォータニオンの導関数
    def qua_dot(self, omega, q):
        omega = self.quaternion(omega)        # omega:機軸座標系の角速度ベクトル
        return -0.5 * self.product(omega, q)
