ARCH?=cpu
MACHINE?=my_laptop
CODE_NAME?=phi4_2d

ifeq ($(MACHINE),studs)
#{


ifeq ($(ARCH),cpu)
#{
BACKEND=kernel_cpu.o
CC=cc
CFLAGS=--std=c11 -O3 -fopenmp -D ARCH=0 -D MACHINE=0
EXEC=${CODE_NAME}_cpu
#}

else

ifeq ($(ARCH),mpi)

#{
BACKEND=kernel_mpi.o comm_mpi.o
CC=/usr/lib/hpc/gnu7/mpi/openmpi/3.1.6/bin/mpicc
CFLAGS=-fopenmp -D ARCH=1 -D MACHINE=0 -O3 --std=c11
EXEC=${CODE_NAME}_mpi
#}

endif
endif

CC_F2I=cc
CFLAGS_F2I=--std=c11 -O3 -fopenmp -D ARCH=0 -D MACHINE=0

#}
#-----------------------------
else

ifeq ($(MACHINE),compass)
#{

ifeq ($(ARCH),cpu)
#{
BACKEND=kernel_cpu.o
CC=cc
CFLAGS=--std=c11 -O3 -fopenmp -D ARCH=0 -D MACHINE=0 
#CFLAGS=-ggdb
EXEC=${CODE_NAME}_cpu
#}
else
#{
ifeq ($(ARCH),mpi)
BACKEND=kernel_mpi.o comm_mpi.o
CC=mpicc
CFLAGS=-fopenmp -D ARCH=1 -D MACHINE=0 -O3 --std=c11 
EXEC=${CODE_NAME}_mpi
endif
#}
endif
CC_F2I=cc
CFLAGS_F2I=--std=c11 -O3 -fopenmp -D ARCH=0 -D MACHINE=0 -I/usr/local/anaconda3/include/libpng16 -L/usr/local/anaconda3/lib/

#}
#-----------------------------

else

ifeq ($(MACHINE),my_laptop)
#{

ifeq ($(ARCH),cpu)
#{
BACKEND=kernel_cpu.o
CC=cc
CFLAGS=--std=c11 -O3 -fopenmp -D ARCH=0 -D MACHINE=0
#CFLAGS=-ggdb
EXEC=${CODE_NAME}_cpu
#}
else
#{
ifeq ($(ARCH),mpi)
BACKEND=kernel_mpi.o comm_mpi.o
CC=mpicc
CFLAGS=-fopenmp -D ARCH=1 -D MACHINE=0  -O3 --std=c11 
EXEC=${CODE_NAME}_mpi
endif
#}
endif

CC_F2I=cc
CFLAGS_F2I=--std=c11 -O3 -fopenmp -D ARCH=0 -D MACHINE=0

#}
endif
endif
endif

HOSTOBJECTS=physics.o vanilla.o random.o lecuyer.o kernel_global.o  timer.o
HEADERS=kernel.h vanilla.h random.h lecuyer.h macros.h timer.h comm.h

.PHONY: all clean executable

all: executable field2image
#all: executable 

executable: $(EXEC)

clean:
	rm -f $(HOSTOBJECTS) $(BACKEND) $(EXEC) field2image

$(EXEC): $(HOSTOBJECTS) $(BACKEND)
	$(CC) $(CFLAGS) $(HOSTOBJECTS) $(BACKEND) $(LDFLAGS) -lm -o $@

field2image: field2image.c
	$(CC_F2I) $(CFLAGS_F2I) $< -lpng -lm -o $@

kernel_cpu.o: kernel_cpu.c macros.h kernel.h
	$(CC) $(CFLAGS) -c $< -o $@

kernel_gpu.o: kernel_gpu.cu macros.h kernel.h   
	$(NVCC) $(CUDA_CFLAGS) $(CUDA_INCFLAGS) $(LDFLAGS) -c $< -o $@

kernel_mpi.o:kernel_mpi.c macros.h kernel.h
	$(CC) $(CFLAGS) -c $< -o $@

comm_mpi.o:comm_mpi.c comm.h 
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@


