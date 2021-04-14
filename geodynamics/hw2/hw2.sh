#!/bin/sh

gmt begin elevation pdf
gmt xyz2grd elevation.xyz -Rg -I1/1 -Gelevation.grd -r -V
gmt grdimage elevation.grd -Rd -JN6i -Bx60 -By30 -Cgeo
gmt colorbar -DJBC+w10c/0.25c+o0.0c/1c -Bxa3f+l"Elevation/km"
gmt end show

gmt begin diffence pdf
gmt xyz2grd diffence.xyz -Rg -I1/1 -Gdiffence.grd -r -V
gmt grdimage diffence.grd -Rd -JN6i -Bx60 -By30 -Cgeo
gmt colorbar -DJBC+w10c/0.25c+o0.0c/1c -Bxa3f+l"Difference/km"
gmt end show

# local 

gmt begin e_local pdf
gmt xyz2grd elevation.xyz -Rg -I1/1 -Gelevation.grd -r -V
gmt grdimage elevation.grd -R-80/-30/-60/10 -JM15c -Bx60 -By30 -Cgeo
gmt colorbar -DJBC+w10c/0.25c+o0.0c/1c -Bxa3f+l"Elevation/km"
gmt end show

gmt begin d_local pdf
gmt xyz2grd diffence.xyz -Rg -I1/1 -Gdiffence.grd -r -V
gmt grdimage diffence.grd -R-80/-30/-60/10 -JM15c -Bx60 -By30 -Cgeo
gmt colorbar -DJBC+w10c/0.25c+o0.0c/1c -Bxa3f+l"Difference/km"
gmt end show
