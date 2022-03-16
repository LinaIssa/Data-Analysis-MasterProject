#!/bin/bash
# -*- coding: utf-8 -*-

#
# 2017.03.03 @dlan

STATIC=0
if [ "$1" == "-s" ]; then
    echo "enabling static mode"
    STATIC=1
    shift;
fi

CONFIN=$1

switch_static(){
   cp ${CONFIN} ${CONFIN}_init
   
   # Mac (bsd) version
   # sed -E -e 's/([^ ]+)(.+specfile.+)/\1_init\2/' -i ''  ${CONFIN}_init
   # sed -E -e 's/([^ ]+)(.+chemfile.+)/\1_init\2/' -i ''  ${CONFIN}_init

   # sed -E -e 's/[CSJ](.+shock_type.+)/S\1/' -i ''  ${CONFIN}_init
   # sed -E -e 's/[123](.+Nfluids.+)/1\1/' -i ''  ${CONFIN}_init
   # sed -E -e 's/[^ ]+(.+timeJ.+)/1.00E+09\1/' -i '' ${CONFIN}_init
   # sed -E -e 's/[^ ]+(.+duration_max.+)/1.00E+09\1/' -i '' ${CONFIN}_init
   
   # Linux (gnu) version
   sed -r -e 's/([^ ]+)(.+specfile.+)/\1_init\2/' -i  ${CONFIN}_init
   sed -r -e 's/([^ ]+)(.+chemfile.+)/\1_init\2/' -i  ${CONFIN}_init

   sed -r -e 's/[CSJ](.+shock_type.+)/S\1/' -i  ${CONFIN}_init
   sed -r -e 's/[123](.+Nfluids.+)/1\1/'  -i  ${CONFIN}_init
   sed -r -e 's/[^ ]+(.+timeJ.+)/1.00E+09\1/' -i  ${CONFIN}_init
   sed -r -e 's/[^ ]+(.+duration_max.+)/1.00E+09\1/' -i  ${CONFIN}_init
   
   cp ${CONFIN}_init ${CONFIN}
}

if [ "$STATIC" -eq "1" ]; then
    echo "building static .in..."
    echo " backuping original .in file"
    cp ${CONFIN} ${CONFIN}.orig
    echo " switch static"
    switch_static
    echo "run mhd_vode with the static .in file..."
    ./mhd_vode ${CONFIN}
    echo "restore original .in file..."
    cp ${CONFIN}.orig ${CONFIN}
else
    echo "run mhd_vode..."
   ./mhd_vode ${CONFIN}
fi

exit 0
