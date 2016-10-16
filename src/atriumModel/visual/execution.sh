#!/bin/bash

# Script to permits the recolection of time information, the first parameter
# is the name of the file containing the time
# the second parameter is the number of executions to do
# the third is the number of cells in x dimension
# the fourth is the number of cells in y dimension
fileName=$1
numberExecutions=$2
xCellNumber=$3
yCellNumber=$4

touch ${fileName}
echo 'user,system,total' > ${fileName}

for ((i=0;i<${numberExecutions};i++))
do
    /usr/bin/time -f'%U,%S,%e' -o ${fileName} -a ./build/visual ${xCellNumber} ${yCellNumber} 0 ../pythonScripts/nuevoTest.py
    echo $i
done
