all: lib serial_FEM
	- ./serial_FEM
	
lib:
	( cd ../../Lib ; $(MAKE) )
	
serial_FEM:
	gcc -I../../Include -g serial_FEM ../../Lib/libphc.a -o serial_FEM
