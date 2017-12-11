# -*- coding: utf-8 -*-


class Wind(object):
    def __init__(self, z_R, v_R):
        self.z_R = z_R
        self.v_R = v_R

    def wind(self, hight):
        # hight<0のとき
        hight = max(0, hight)
        return self.v_R * ((hight / self.z_R)**(1/7.0))
