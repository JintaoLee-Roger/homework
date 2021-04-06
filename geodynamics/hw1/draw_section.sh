#!/bin/sh
gmt set FONT_ANNOT_PRIMARY=10
gmt set FONT_LABEL=12

gmt begin section pdf
gmt basemap -R0/12/-8000/2000 -JX10c/3c -BWrtb -Bya2000f+l"Elevation (m)"
# 生成地震深度颜色表
echo 0 purple@30 70 purple@30 > depth.cpt
echo 70 green@30 300 green@30 >> depth.cpt
echo 300 red@30 700 red@30 >> depth.cpt

# 标注 A B 位置
echo 0 2200 A > tmp
echo 12 2200 B >> tmp
gmt text tmp -F+f10p+jBC -N -D0c/0.1c

# 地形投影
gmt project -C117/6 -E129/6 -G0.1 | gmt grdtrack -G@earth_relief_03m > tmp
echo 0 0 > tmp2 
echo 12 0 >> tmp2
gmt plot tmp2 -Wblack -Glightblue -L+y-8000
gmt plot tmp -i2,3 -Wblack -Ggray -L+y-8000

# benioff zone
gmt basemap -R0/12/0/700 -JX10c/-5c -Bya200f100+l"Focal depth (km)" -Bxa2f1+l"Distance"+u"\260" -BWSrt -Y-5.5c
echo 5.5 480 Benioff zone | gmt text -F+f14p,5,blue=solid+a50+jBL

gmt coupe phili -Q -L -Sm0.3c -Aa117/6/129/6/90/300/0/700f -Zdepth.cpt

( 
cat << EOF 
H 8,4 Events (Mw >= 5)
S 0.1i c 0.20 purple 0.1p,black 0.20i 0-70 km
S 0.1i c 0.20 green 0.1p,black 0.20i 70-300 km
S 0.1i c 0.20 red 0.1p,black 0.20i 300-700 km
EOF
) > tmp

gmt legend tmp -DjBR+w1i+l1.5+o0.05i/0.04i -F+g255+p0.25p

gmt end show

rm *.cpt tmp*