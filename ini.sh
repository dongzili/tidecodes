#!/bin/bash

sed -i 's/kperp0.13/l300/g' *.cpp
sed -i 's/\[80\]/\[\]/g' *.cpp
sed -i 's/kc0.5/kc0.6/g' *.cpp
