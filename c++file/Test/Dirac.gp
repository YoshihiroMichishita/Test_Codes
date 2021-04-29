L="w l lw 0.3"

unset tics 
unset border
unset key

set xrange[-0.2:0.2]
set yrange[-0.2:0.2]

set ticslevel 0
set view 68, 30, 1, 1

set out "Dirac.png" # ファイル名を指定
set terminal pngcairo size 600, 600 # アンチエイリアスpng形式 (600*600 px) で出力
splot 'Dirac.dat' u 1:2:3 @L lc 7,'Dirac.dat' u 1:2:4 @L lc 6