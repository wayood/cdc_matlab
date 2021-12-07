# Purpose
middle planでの軌道補正アルゴリズムを適用させる。
# Usage
cdc_middle_save.m -> global(Astar)でのパスを生成、重いためここで分割  
middle_planner_ex.m -> パスを用いてDWAで走行、ここで軌道を保存  
middle_planner_cdc_ex.m -> 軌道補正入れたバージョン、あとは上行と同様  
curvature.m,circument.m -> 曲率計算、上二行で使用  
plot_potential_cdc.m -> ポテンシャル場で安全性を評価  
Potential_Field.m -> ポテンシャル場を生成  
# Memo
評価基準あいまい