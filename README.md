# Purpose
軌道補正の事前に破綻検知をできること,別の問題を発見する
# Usage
```
>> setup
```
それぞれの使用用途　
```
cdc_middle_save.m -> global(Astar)でのパスを生成、重いためここで分割  
middle_planner_ex.m -> パスを用いてDWAで走行、ここで軌道を保存  
middle_planner_cdc_ex.m -> 軌道補正入れたバージョン、あとは上行と同様  
DouglasPeuker.m -> 折れ線近似による補正ポイントの制御  
plot_potential_cdc.m -> ポテンシャル場で安全性を評価  
Potential_Field.m -> ポテンシャル場を生成  
cdc_teleoperation.m -> WPを人が指定し,その間を機械計算で走行する.並列処理により経路追加しながら走行ができる状態
cdc_teleoperation_waypoint_add.m -> ボロノイ図で破綻検知後やローカルミニマムに対して解決中
cdc_teleoperation_autonomous.m -> 破綻検知後のボロノイ図を自律軌道計画として導入  
```
# Plan proccessing　
解析は出来ている前提で進む.それ以降のCDCによる長距離,長寿命を考える.
