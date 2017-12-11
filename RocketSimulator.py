# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import g
from scipy.integrate import odeint
import pandas as pd
import pyproj
import simplekml
from quaternion import Quaternion
from mass import Mass
from force import Force
from wind import Wind
import angle_of_attack as aoa


class RocketSimulator(object):

    def __init__(self):
        self.T = 150
        self.dt = 0.01
        self.time = np.arange(0.0, self.T, self.dt)
        self.N = len(self.time)

    def initialize(self, design, launch_condition):
        """ 初期化 """
        self.name = design['name']
        self.m_af = design['m_af']
        self.I_af = design['I_af']
        self.CP = design['CP']
        self.CG_a = design['CG_a']
        self.d = design['d']
        self.area = np.pi * (self.d ** 2) / 4.0
        self.len_a = design['len_a']
        self.inertia_z0 = design['inertia_z0']
        self.inertia_zT = design['inertia_zT']
        self.engine = design['engine']
        self.me_total = design['me_total']
        self.me_prop = design['me_prop']
        self.len_e = design['len_e']
        self.d_e = design['d_e']

        self.p0 = np.array([0., 0., 0.])  # position(x, y, z)
        self.condition_name = launch_condition['name']
        self.theta0 = launch_condition['AngleOfFire']
        self.phi0 = launch_condition['azimuthal']
        self.launch_rod = launch_condition['launch_rod']
        self.v0 = np.array([0., 0., 0.])  # velocity(vx, vy, vz)
        self.ome0 = np.array([0., 0., 0.])
        self.density = launch_condition['density']
        self.wind_R = launch_condition['StandardWind']
        self.z_R = launch_condition['StandardHeight']
        self.beta = launch_condition['WindDirection']  # wind direction
        self.wind_direction = np.array([np.cos(self.beta), np.sin(self.beta), 0.0])
        self.qua_theta0 = np.array([np.cos(0.5 * self.theta0), np.sin(0.5 * self.theta0), 0., 0.])  # x軸theta[rad]回転, 射角
        self.qua_phi0 = np.array([np.cos(0.5 * self.phi0), 0., 0., np.sin(0.5 * self.phi0)])  # z軸phi[rad]回転, 方位角
        self.wind_direction = np.array([np.cos(self.beta), np.sin(self.beta), 0.0])

        self.engine_data = np.loadtxt(self.engine)

        self.force = Force(self.area, self.engine_data, self.T, self.density)
        self.thrust = self.force.thrust()

        self.mass = Mass(self.m_af, self.I_af, self.CG_a, self.len_a, self.inertia_z0, self.inertia_zT, self.me_total,
                         self.me_prop, self.len_e, self.d_e, self.force.burn_time, self.T)
        self.M = self.mass.mass()
        self.Me = self.mass.me_t()
        self.Me_dot = self.mass.me_dot()
        self.CG = self.mass.CG()
        self.CG_dot = self.mass.CG_dot()
        self.Ie = self.mass.iexg()
        self.Inertia = self.mass.inertia()
        self.Inertia_z = self.mass.inertia_z()
        self.Inertia_dot = self.mass.inertia_dot()
        self.Inertia_z_dot = self.mass.inertia_z_dot()

        self.wind = Wind(self.z_R, self.wind_R)

    def deriv(self, pi, vi, quai, omei, t):
        """ 運動方程式 """
        qt = Quaternion()
        # 機軸座標系の推力方向ベクトル
        r_Ta = np.array([0., 0., 1.0])
        # 慣性座標系重力加速度
        gra = np.array([0., 0., -g])
        # 機軸座標系の空力中心位置
        r = np.array([0., 0., self.CG(t) - self.CP])
        # 慣性座標系の推力方向ベクトル
        r_T = qt.rotation(r_Ta, qt.coquat(quai))
        r_T /= np.linalg.norm(r_T)
        # 慣性テンソル
        I = np.diag([self.Inertia(t), self.Inertia(t), self.Inertia_z(t)])
        # 慣性テンソルの時間微分
        I_dot = np.diag([self.Inertia_dot(t), self.Inertia_dot(t), self.Inertia_z_dot(t)])
        # 慣性座標系対気速度
        v_air = self.wind.wind(pi[2]) * self.wind_direction - vi
        # 迎角
        alpha = aoa.aoa(qt.rotation(v_air, quai))
        # ランチロッド垂直抗力
        N = 0
        # ランチロッド進行中
        if np.linalg.norm(pi) <= self.launch_rod and r_T[2] >= 0:
            Mg_ = self.M(t) * gra - np.dot(self.M(t) * gra, r_T) * r_T
            D_ = self.force.drag(alpha, v_air) - np.dot(self.force.drag(alpha, v_air), r_T) * r_T
            N = -Mg_ - D_
        # 慣性座標系加速度
        v_dot = gra + (self.thrust(t) * r_T + self.force.drag(alpha, v_air) + N) / self.M(t)
        # クォータニオンの導関数
        qua_dot = qt.qua_dot(omei, quai)
        # 機軸座標系角加速度
        ome_dot = np.linalg.solve(I, - np.cross(r, qt.rotation(self.force.drag(alpha, v_air), quai))
                                  - np.dot(I_dot, omei) - np.cross(omei, np.dot(I, omei)))
        # ランチロッド進行中
        if np.linalg.norm(pi) <= self.launch_rod:
            # ランチロッド進行中は姿勢が一定なので角加速度0とする
            ome_dot = np.array([0., 0., 0.])

        return vi, v_dot, qua_dot, ome_dot

    def simulate(self, method='RungeKutta', log=False):
        """ 数値計算 """
        qt = Quaternion()
        p = np.empty((self.N+1, 3))
        v = np.empty((self.N+1, 3))
        qua = np.empty((self.N+1, 4))
        ome = np.empty((self.N+1, 3))
        p[0] = self.p0
        v[0] = self.v0
        qua[0] = qt.product(self.qua_phi0, self.qua_theta0)
        ome[0] = self.ome0
        count = 0

        for (i, t) in enumerate(self.time):
            if method == 'RungeKutta':
                # Runge-Kutta method
                pk1, vk1, quak1, omek1 = self.deriv(p[i], v[i], qua[i], ome[i], t)
                pk2, vk2, quak2, omek2 = self.deriv(p[i] + pk1 * self.dt * 0.5, v[i] + vk1 * self.dt * 0.5, qua[i] + quak1 * self.dt * 0.5, ome[i] + omek1 * self.dt * 0.5, t + self.dt * 0.5)
                pk3, vk3, quak3, omek3 = self.deriv(p[i] + pk2 * self.dt * 0.5, v[i] + vk2 * self.dt * 0.5, qua[i] + quak2 * self.dt * 0.5, ome[i] + omek2 * self.dt * 0.5, t + self.dt * 0.5)
                pk4, vk4, quak4, omek4 = self.deriv(p[i] + pk3 * self.dt, v[i] + vk3 * self.dt, qua[i] + quak3 * self.dt, ome[i] + omek3 * self.dt, t + self.dt)
                p[i + 1] = p[i] + (pk1 + 2 * pk2 + 2 * pk3 + pk4) * self.dt / 6
                v[i + 1] = v[i] + (vk1 + 2 * vk2 + 2 * vk3 + vk4) * self.dt / 6
                qua[i + 1] = qua[i] + (quak1 + 2 * quak2 + 2 * quak3 + quak4) * self.dt / 6
                ome[i + 1] = ome[i] + (omek1 + 2 * omek2 + 2 * omek3 + omek4) * self.dt / 6

            elif method == 'Euler':
                # Euler method
                p_dot, v_dot, qua_dot, ome_dot = self.deriv(p[i], v[i], qua[i], ome[i], t)
                p[i + 1] = p[i] + p_dot * self.dt
                v[i + 1] = v[i] + v_dot * self.dt
                qua[i + 1] = qua[i] + qua_dot * self.dt
                ome[i + 1] = ome[i] + ome_dot * self.dt

            if t <= self.force.burn_time:
                p[i + 1][2] = max(0., p[i + 1][2])

            # vz<0かつz<0のとき計算を中断
            if v[i + 1][2] < 0 and p[i + 1][2] < 0:
                count = i + 1
                print("Calculation was successfully completed.")
                break

        self.p = p[:count + 1]
        self.v = v[:count + 1]
        self.qua = qua[:count + 1]
        self.ome = ome[:count + 1]
        if log:
            np.savetxt(self.name + "_" + self.condition_name + "_position.csv", p, delimiter=",")
            print("Position file was created.")

        print("Altitude is " + str(max(p[:, 2])) + "[m].")

    def output_kml(self, place):
        # 原点からの距離
        def dis2d(x, y):
            return pow(pow(x, 2) + pow(y, 2), 0.5)

        # y軸とベクトル(x, y)のなす角
        def ang2d(x, y):
            # y軸上
            if x == 0:
                if y >= 0:
                    return 0.0
                else:
                    return 180.0
            # x軸上
            if y == 0:
                if x >= 0:
                    return 90.0
                else:
                    return 270.0
            # 第1象限
            if x > 0 and y > 0:
                return np.arctan(x / y) * 180 / np.pi
            # 第2象限
            if x > 0 and y < 0:
                return 180.0 + np.arctan(x / y) * 180 / np.pi
            # 第3象限
            if x < 0 and y < 0:
                return 180.0 + np.arctan(x / y) * 180 / np.pi
            # 第4象限
            if x < 0 and y > 0:
                return 360.0 + np.arctan(x / y) * 180 / np.pi

        distance = [dis2d(self.p[i, 0], self.p[i, 1]) for i in range(len(self.p))]

        angle = [ang2d(self.p[i, 0], self.p[i, 1]) for i in range(len(self.p))]

        coordinate0 = place[1]
        latitude = coordinate0[0]
        longitude = coordinate0[1]
        geod = pyproj.Geod(ellps='WGS84')
        newLong = []
        newLat = []
        invAngle = []
        for i in range(len(self.p)):
            nlon, nlat, nang = geod.fwd(longitude, latitude, angle[i], distance[i])
            newLong.append(nlon)
            newLat.append(nlat)
            invAngle.append(nang)

        kml = simplekml.Kml(open=1)
        cood = []
        for i in range(len(self.p)):
            cood.append((newLong[i], newLat[i], self.p[i, 2]))
        ls = kml.newlinestring(name=self.name+"'s Path")
        ls.coords = cood
        ls.style.linestyle.width = 3
        ls.style.linestyle.color = simplekml.Color.blue
        ls.altitudemode = simplekml.AltitudeMode.relativetoground
        ls.extrude = 0
        kml.save(self.name + "_" + place[0] + '_' + self.condition_name + "_path.kml")
        print("Kml file for Google Earth was created.")


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    design = pd.read_csv('design.csv')
    condition = pd.read_csv('condition.csv')

    rs = RocketSimulator()
    rs.initialize(design.loc[0], condition.loc[0])
    rs.simulate(method='RungeKutta', log=True)
    
    places = [["NOSHIRO_A", [40.138159, 139.984342]],
              ["Rits_Ground1", [34.977867, 135.963783]],
              ["Kashiwa", [35.897433, 139.937890]]]
    rs.output_kml(places[2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(rs.p[:, 0], rs.p[:, 1], rs.p[:, 2])
    ax.set_xlabel(u"x")
    ax.set_ylabel(u"y")
    ax.set_zlabel(u"z")
    ax.set_zlim([0., 250])
    plt.show()

    # plt.figure(figsize=(6, 8))
    # plt.plot(r2s.p[:, 1], r2s.p[:, 2])
    # plt.ylim(0, 250)
    # plt.xlabel('y')
    # plt.ylabel('z')
    # plt.show()
