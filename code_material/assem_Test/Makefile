all: lib assem_A
	- gdb ./assem_A
	
run: lib assem_A
	- ./assem_A
	
lib:
	( cd ../Lib ; $(MAKE) )
	
assem_A:
	gcc -I../Include -g assem_A.c ../Lib/libhpc.a -o assem_A -lm
