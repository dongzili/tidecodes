
#g++ -g -Wall tide.cpp smooth.cpp shear.cpp k3d_noisy.cpp powerspec.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -o tide00
#g++ -g -Wall plotsmooth.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -o tide00

#g++ -g -Wall smooth.cpp -L/opt/fftw-3.3.3/lib -I/opt/fftw-3.3.3/include -lfftw3 -lm
nohup ./tides03 > log_tides03 2>&1 &
sleep 3600
nohup ./tides04 > log_tides04 2>&1 &
sleep 3600
nohup ./tides05 > log_tides05 2>&1 &
#vfield.cpp filterv.cpp momenfield.cpp powerv.cpp
