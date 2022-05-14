
CC = gcc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib hpc_demo1 hpc_demo2 
	- ./hpc_demo1 problem1 2
	- ./hpc_demo1 problem1 4


lib:
	( cd ../Lib ; $(MAKE) )

hpc_demo1: lib hpc_demo1.c Makefile
	$(CC) $(CF) $(I) -o hpc_demo1 hpc_demo1.c $(CS)

hpc_demo2: lib hpc_demo2.c Makefile
	$(CC) $(CF) $(I) -o hpc_demo2 hpc_demo2.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo1 hpc_demo2 hpc_demo3 hpc_demo4 hpc_demo5 \
    *.a *.dSYM *.obj *.dll
