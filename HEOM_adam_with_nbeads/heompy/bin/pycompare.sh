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
dir=${prog}_${timestamp}_${name}
mkdir -p ${dir}
cd ${dir}
cp ${inputdir}/${input} .
heompy $input $output 2>&1 | tee $stdout &
heom ${input} | tee ${stdout}_fort &
wait
read -p "Do you want to plot the graph?" yn
case $yn in
    [Yy]* ) plot.py $output out.dat ;;
    [Nn]* ) ;;
esac
