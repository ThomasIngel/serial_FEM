all: lib parallel
	
lib:
	( cd ../code_material/Lib ; $(MAKE) )
	
parallel:
	mpicc -I../code_material/Include -g parallel.c ../code_material/Lib/libhpc.a -lm -o parallel
