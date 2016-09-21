
#g++ -g -Wall tide.cpp smooth.cpp shear.cpp k3d_noisy.cpp powerspec.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -o tide00
#g++ -g -Wall plotsmooth.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -o tide00

#g++ -g -Wall smooth.cpp -L/opt/fftw-3.3.3/lib -I/opt/fftw-3.3.3/include -lfftw3 -lm
nohup ./newtides03 > log_newtides03 2>&1 &
nohup ./newtides04 > log_newtides04 2>&1 &
nohup ./newtides05 > log_newtides05 2>&1 &
#vfield.cpp filterv.cpp momenfield.cpp powerv.cpp
