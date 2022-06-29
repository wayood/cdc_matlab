# Purpose
middle planでの軌道補正アルゴリズムを適用させる。
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
cdc_teleoperation -> WPを人が指定し,その間を機械計算で走行する.  
```
# Plan proccessing
評価の距離トレードオフを議論　
path planningとteleoperationとの統合->シームレス化を促す　
解析は出来ている前提で進む.それ以降のCDCによる長距離,長寿命を考える.
