F90 =  gfortran -fbounds-check 

all  :  obj global 
	
global: 
	$(F90) -o run b-Bayesian_binned.f90 \
	-L./ -lroutines -lm
obj:
	$(F90) -c src_likelihood/sort_functions.f90
	$(F90) -c src_likelihood/nrtype.f90
	$(F90) -c src_likelihood/nrutil.f90
	$(F90) -c src_likelihood/nr.f90
	$(F90) -c src_likelihood/likelihoods_FMD_gridsearch.f90
	$(F90) -c src_likelihood/likelihoods_FMD_random_complete.f90

	ar -r libroutines.a sort_functions.o nrtype.o nrutil.o nr.o likelihoods_FMD_gridsearch.o likelihoods_FMD_random_complete.o
	\rm ./*.o
clean:
	/bin/rm *.a
	
