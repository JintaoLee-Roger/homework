#!/bin/sh
gmt set FONT_ANNOT_PRIMARY=10
gmt set FONT_LABEL=14

# gmt begin globe_relief pdf
# gmt grdimage ETOPO1_Bed_g_gdal.grd -JH180/15c -I+d -Cgeo
# gmt colorbar -DJBC+w10c/0.25c+o0.0c/0.3c -Bxa3f+l"Topography/km" -C+Ukm
# gmt end show

gmt begin globe_relief pdf
gmt grdimage @earth_relief_03m -JH180/15c -I+d
gmt colorbar -DJBC+w10c/0.25c+o0.0c/0.3c -Bxa3f+l"Topography/km" -C+Ukm
gmt end show

gmt begin phlippines_relif pdf
gmt grdimage @earth_relief_15s -R115/130/0/22 -JM15c -I+d -Ba -BNWes
gmt colorbar -DJBC+w10c/0.25c+o0.0c/0.3c -Bxa3f+l"Topography/km" -C+Ukm
gmt end show