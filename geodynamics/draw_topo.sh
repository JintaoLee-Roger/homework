#!/bin/sh

gmt begin global_relief pdf
gmt grdimage @earth_relief_03m -JH180/10c
gmt end show
