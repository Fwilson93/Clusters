module load cray-fftw
gfortran test.f95 -I/home/n03/n03/fwilson/SHTOOLS/include -fopenmp -m64 -fPIC -O3 -std=gnu -ffast-math -L/home/n03/n03/fwilson/SHTOOLS/lib -lSHTOOLS-mp -o clstr
