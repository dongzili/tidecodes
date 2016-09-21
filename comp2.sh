
#g++ -g -Wall tide.cpp smooth.cpp shear.cpp k3d_noisy.cpp powerspec.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -o tide00
#g++ -g -Wall plotsmooth.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -o tide00

g++ -g -Wall powerspec.cpp momenorigin.cpp filterv.cpp powerv.cpp momenfield.cpp powermomen.cpp filter_neg.cpp tide.cpp smooth_fs.cpp recon.cpp vfield.cpp -L /opt/fftw-3.3.3/lib -I /opt/fftw-3.3.3/include -lfftw3 -lm -L/opt/gsl-1.15-gcc-4.4.6/lib -I/opt/gsl-1.15-gcc-4.4.6/include -lgsl -lgslcblas -o newtides00


#g++ -g -Wall smooth.cpp -L/opt/fftw-3.3.3/lib -I/opt/fftw-3.3.3/include -lfftw3 -lm
#vfield.cpp filterv.cpp momenfield.cpp powerv.cpp
