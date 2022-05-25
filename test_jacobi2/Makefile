all: lib test_jacobi
	- ./test_jacobi
	
lib:
	( cd ../code_material/Lib ; $(MAKE) )
	
test_jacobi:
	gcc -I../code_material/Include -g test_jacobi_seriell_v2.c ../code_material/Lib/libhpc.a -lm -o test_jacobi
