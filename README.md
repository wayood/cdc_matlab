# Purpose
middle planでの軌道補正を適用し、軌道補正自体と評価全体の問題点を見つける。
# Usage
cdc_middle_save.m -> global(Astar)でのパスを生成、重いためここで分割 
cdc_middle_dwa_save.m ->global(DWA)で引き続けるもの　
middle_planner_ex.m -> 階層パスプランニング  
middle_planner_cdc_ex.m -> 軌道補正入れたバージョン、あとは上行と同様  
Douglas Peucker.m -> 折線近似によるWPの自動化　
plot_potential_cdc.m -> ポテンシャル場で安全性を評価  
Potential_Field.m -> ポテンシャル場を生成  
# Memo
現状　
50,100,150,200で走行を確認　
評価について　
手法の改善について考えるべき、評価を見てから姿勢情報を追加すること