#!/bin/sh
gmt set FONT_ANNOT_PRIMARY=10
gmt set FONT_LABEL=12

gmt begin earthquake_globe pdf
gmt grdimage @earth_relief_03m -I+d -Ba -BNWes
# 生成地震深度颜色表
echo 0 purple@30 70 purple@30 > depth.cpt
echo 70 green@30 300 green@30 >> depth.cpt
echo 300 red@30 800 red@30 >> depth.cpt
gmt meca cmt -Sm0.3c -Zdepth.cpt #-JH180/15c
# 画图例 下方
( 
cat << EOF 
H 10,0 Earthquake Location
L 6,0 C (2011/01/01-2021/01/01)
G 0.1c
L 8,0 C Magnitude(Mw)
G 0.1c
N 3
S 0.1i c 0.11i - black 0.2i 6
S 0.1i c 0.12i - black 0.2i 7
S 0.1i c 0.13i - black 0.2i 8
N 1
G 0.1c
L 8,0 C Depth(km)
G 0.1c
N 3
L 7,0 R @;purple;0-70km@;;
L 7,0 R @;green;70-300km@;;
L 7,0 R @;red;300-800km@;;
G 1.0c
B globe 0.3i 0.08i+ml -Bxa3f+l"Topo(km)" -C+Ukm --FONT_ANNOT_PRIMARY=8p --MAP_FRAME_WIDTH=1p --FONT_LABEL=10p
EOF
) > tmp
gmt legend tmp -F+gazure1@10 -DjBC+w8c+l1.2+o0.0c/-5c -C0.1i/0.1i
echo 123 33.05 43 3.62 -0.44 -3.18 0.90 2.46 -1.35 24 0 0 > tmp
echo 126.6 33.05 171 -0.71 -0.26 0.96 0.44 0.81 -0.07 24 0 0 >> tmp
echo 130.5 33.05 302 0.34 0.16 -0.50 -0.77 -4.57 -1.58 24 0 0 >> tmp    
gmt meca tmp -Sm0.3c -Zdepth.cpt
gmt end show

gmt begin earthquake_phlippins pdf
# 画地形起伏图
gmt grdimage @earth_relief_15s -R115/130/0/22 -JM15c -I+d -Ba -BNWes
# 生成地震深度颜色表
echo 0 purple@30 70 purple@30 > depth.cpt
echo 70 green@30 300 green@30 >> depth.cpt
echo 300 red@30 800 red@30 >> depth.cpt
# 画beach ball
gmt meca phili -Sm0.3c -Zdepth.cpt -JM15c 
# 画侧线
echo 117 6 A > tmp
echo 129 6 B >> tmp
gmt plot tmp -W1p,black,-.-
gmt text tmp -F+f10p -D0c/-0.2c 
# 画图例 (右上角)
( 
cat << EOF 
H 10,0 Earthquake Location
L 6,0 C (2011/01/01-2021/01/01)
G 0.1c
L 8,0 C Magnitude(Mw)
G 0.1c
N 4
S 0.1i c 0.10i - black 0.2i 5
S 0.1i c 0.11i - black 0.2i 6
S 0.1i c 0.12i - black 0.2i 7
S 0.1i c 0.13i - black 0.2i 8
N 1
G 0.1c
L 8,0 C Depth(km)
G 0.1c
N 3
L 7,0 R @;purple;0-70km@;;
L 7,0 R @;green;70-300km@;;
L 7,0 R @;red;300-800km@;;
N 1
G 0.1c
L 8,0 C Survey Line
G 0.1c
N 2
L 7,0 R A(117 6)
L 7,0 R B(129 6)
EOF
) > tmp
gmt legend tmp -F+gazure1@10 -DjTR+w4c+l1.2 -C0.1i/0.1i
echo 123 33.05 43 3.62 -0.44 -3.18 0.90 2.46 -1.35 24 0 0 > tmp
echo 126.6 33.05 171 -0.71 -0.26 0.96 0.44 0.81 -0.07 24 0 0 >> tmp
echo 130.5 33.05 302 0.34 0.16 -0.50 -0.77 -4.57 -1.58 24 0 0 >> tmp    
gmt meca tmp -Sm0.3c -Zdepth.cpt
gmt colorbar -DJBC+w10c/0.25c+o0.0c/0.3c -Bxa3f+l"Topography/km" -C+Ukm
gmt end show