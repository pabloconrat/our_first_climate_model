gcc -Wall -c test.c -I .
gcc -Wall -c ascii.c -I .
g++ -Wall -std=gnu++14 -c repwvl_thermal.cc
g++ -o test test.o repwvl_thermal.o ascii.o -lnetcdf_c++4 -L/opt/local/lib
