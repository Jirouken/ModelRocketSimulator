# ModelRocketSimulator
モデルロケットの設計データと打ち上げ条件から軌道を計算する．

## description
ロケットの運動方程式をRunge-Kutta4/Euler法により数値計算を行う．
軌道シミュレーション結果をGoogle Earth上に表示するためのkmlファイルを作成する．
![シミュレーション結果](https://github.com/Jirouken/ModelRocketSimulator/blob/master/VESPER_NOSHIRO_A_condition1_path.png)

### 入力
ロケットの設計に基づく物理量 design.csv
- 機体質量
- 機体重心(CG)
- 圧力中心(CP)
- 直径
- 全長
- 機体x,y軸慣性モーメント
- 燃焼前エンジン含む機体z軸慣性モーメント
- 燃焼後エンジン含む機体z軸慣性モーメント
- 燃焼前エンジン質量
- エンジン燃料質量
- エンジン直径
- エンジン全長
- エンジン推力履歴

打ち上げ条件 condition.csv
- 射角
- 方位角
- ランチロッド長さ
- 空気密度
- 風速基準高さ
- 風速
- 風向

## requirement
- Python3
- numpy
- scipy
- pandas
- pyproj
- simplekml

## usage
~~~
from RocketSimulator import RocketSimulator
import pandas as pd

rs = RocketSimulator()
rs.initialize(design.loc[0], condition.loc[0])
rs.simulate(method='RungeKutta', log=True)
place = [["NOSHIRO_A", [40.138159, 139.984342]]
rs.output_kml(place)
~~~
