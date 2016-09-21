
#g++ -g -Wall tide.cpp smooth.cpp shear.cpp k3d_noisy.cpp powerspec.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -o tide00
#g++ -g -Wall plotsmooth.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -o tide00

g++ -g -Wall power2d.cpp  -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -L/opt/gsl-1.15-gcc-4.4.6/lib -I/opt/gsl-1.15-gcc-4.4.6/include -lgsl -lgslcblas -o pow2d00

nohup ./pow2d00 > log_powv2d 2>&1 &

#g++ -g -Wall smooth.cpp -L/opt/fftw-3.3.3/lib -I/opt/fftw-3.3.3/include -lfftw3 -lm
#vfield.cpp filterv.cpp momenfield.cpp powerv.cpp
