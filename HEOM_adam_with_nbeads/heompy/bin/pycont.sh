#!/bin/bash
# Make script stop at first error
set -e
prog=heompy
name=$1
input=input
output=output
stdout=stdout

export PYTHONUNBUFFERED=TRUE

echo "Running pyrun.sh for program " $prog

heompy $input $output 2>&1 | tee $stdout
read -p "Do you want to plot the graph?" yn
case $yn in
    [Yy]* ) plot.py $output ../bath_dense ;;
    [Nn]* ) ;;
esac
