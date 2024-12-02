#!/bin/bash
for directory in */;do
    name=${directory/"/"/""}
    if [ -f ${directory}rho ]; then
        cp ${directory}rho ./rho_${name}
    fi
done
