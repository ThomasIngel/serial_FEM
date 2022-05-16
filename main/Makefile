all: lib serial_FEM
	- ./serial_FEM
	
lib:
	( cd ../code_material/Lib ; $(MAKE) )
	
serial_FEM:
	gcc -I../code_material/Include -g serial_FEM.c ../code_material/Lib/libhpc.a -o serial_FEM
