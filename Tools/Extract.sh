#!/bin/sh

alias extract='tar -xvzf'
temp='--directory ./.Temp'

mkdir ./.Temp

extract Data/Archives/nonrarevariants_exphunter.tar.gz $temp
extract Data/Archives/c9orf72_exphunter.tar.gz $temp
extract Data/Archives/fmr1_exphunter.tar.gz $temp

