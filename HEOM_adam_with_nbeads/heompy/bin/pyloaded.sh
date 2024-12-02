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
cd ${outputdir}
dir=${timestamp}_${name}
mkdir -p ${dir}
cd ${dir}
cp ${inputdir}/${input} .
#dvr_tcf.py $input 2>&1 | tee ${stdout}_dvr
cp ../rho0 .
if [ -f rho0 ]; then
    heompy $input 2>&1 | tee $stdout
fi
read -p "Do you want to plot the graph?" yn
case $yn in
    [Yy]* )
        if [ -f rho ]; then
            plot_rho.py rho -animate
        fi
        if [ -f tcf_qq ]; then
            plot_tcf.py tcf_qq
        fi
        ;;
    [Nn]* ) echo "Not plotting";;
esac
