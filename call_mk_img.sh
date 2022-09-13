#!/usr/bin/bash

imdir=image
mkdir -p ${imdir}
for ((i=0;i<=100000;i+=1000))
do
    is=`printf  %06d  ${i}`
    # ./mk_img.sh     [fname]    [input dir]  [output dir]
    ./mk_img.sh      xyuvp${is}      data        ${imdir}
done
