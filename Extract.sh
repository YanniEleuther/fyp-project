#!/bin/sh

alias extract='tar -xvzf'
temp='--directory ./.Temp'

mkdir ./.Temp

extract Archives/nonrarevariants_exphunter.tar.gz $temp
extract Archives/c9orf72_exphunter.tar.gz $temp
extract Archives/fmr1_exphunter.tar.gz $temp

