all: lib local_sm
	- ./local_sm
	
lib:
	( cd ../code_material/Lib ; $(MAKE) )
	
local_sm:
	gcc -I../code_material/Include -g local_sm.c ../code_material/Lib/libhpc.a -lm -o local_sm
