#!/bin/bash

sed -i 's/\[80\]/\[\]/g' *.cpp
sed -i 's/z1_rpar15_kc0.5/z2_rpar10_kc0.5/g' *.cpp
sed -i 's/\/1\.000/\/2\.000/g' *.cpp
