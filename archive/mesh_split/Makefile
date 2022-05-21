all: lib mesh_split
	- ./mesh_split
	
lib:
	( cd ../code_material/Lib ; $(MAKE) )
	
mesh_split:
	gcc -I../code_material/Include -g mesh_split.c ../code_material/Lib/libhpc.a -lm -o mesh_split
