#!/bin/bash
# Make script stop at first error
set -e
prog=heompy
name=$1
input=input
output=output
stdout=stdout
id=test
outputdir=${OUTDIR}
inputdir=${WORKDIR}/${prog}/input

export PYTHONUNBUFFERED=TRUE

echo "Running pyrun.sh for program " $prog




# Preparing places
timestamp=$(date +%Y-%m-%d_%H-%M-%S)
cd ${outputdir}/outdir
dir=${timestamp}_${name}
dir=${name}
mkdir -p ${dir}
cd ${dir}
cp ${inputdir}/${input} .
touch $timestamp
pypy3 ${WORKDIR}/heompy/bin/heompy $input 2>&1 | tee $stdout
#read -p "Do you want to plot the graph?" yn
#case $yn in
#    [Yy]* )
#        if [ -f tcf* ]; then
#            plot_tcf.py tcf*
#        fi
#        if [ -f rho ]; then
#            plot_rho.py rho -animate
#        fi
#        ;;
#    [Nn]* ) echo "Not plotting";;
#esac
