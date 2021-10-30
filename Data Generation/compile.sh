gfortran -c src/size_dist.f95 -O3
gfortran -c src/miex.f90 -O3
gfortran -c mie_calc.f95 -O3
gfortran -o create_data create_data.f95 mie_calc.f95 src/dmilay_f95.F src/miex.f90 src/size_dist.f95 -O3
