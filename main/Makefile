all: lib serial_FEM
	- mpirun -np 1 ./serial_FEM 1
	
lib:
	( cd ../code_material/Lib ; $(MAKE) )
	
serial_FEM:
	mpicc -I../code_material/Include -g serial_FEM.c ../code_material/Lib/libhpc.a -lm -o serial_FEM
