all: lib test_cg
	- ./test_cg
	
lib:
	( cd ../code_material/Lib ; $(MAKE) )
	
test_cg:
	gcc -I../code_material/Include -g test_cg_seriell_v2.c ../code_material/Lib/libhpc.a -lm -o test_cg
