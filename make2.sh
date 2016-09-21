#!/bin/bash

for k in $( seq 0 5 )
   do

./comp2.sh   
   sleep 4

sed -i 's/\/0'${k}'/\/0'$((k+1))'/g' *.cpp
sed -i 's/tides0'${k}'/tides0'$((k+1))'/g' *.cpp comp2.sh
    sleep 1
done

sed -i 's/\/06/\/00/g' *.cpp
sed -i 's/tides06/tides00/g' *.cpp comp2.sh

