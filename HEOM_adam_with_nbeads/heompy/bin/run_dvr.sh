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

echo "Running pydvr.sh for program dvr_tcf.py"

# Preparing places
timestamp=$(date +%Y-%m-%d_%H-%M-%S)
cd ${outputdir}
dir=${timestamp}_${name}
dir=${name}
mkdir -p ${dir}
cd ${dir}
cp ${inputdir}/${input} .
dvr_tcf.py $input ${name} 2>&1 | tee $stdout
