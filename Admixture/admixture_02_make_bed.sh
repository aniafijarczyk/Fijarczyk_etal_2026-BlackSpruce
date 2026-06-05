#!/bin/bash


ped=admixture_01_input.ped
nameis=$(echo $ped | cut -d"." -f1)
./plink --file ${nameis} --make-bed --out ${nameis} --allow-extra-chr
