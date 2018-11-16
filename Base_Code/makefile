KT= EMPTY
DT= EMPTY
IT= EMPTY

MAKEFLAGS += --jobs=8
CC = nvcc
CCFLAGS = -dc -O4 -D_FORCE_INLINES -DRODC_ON -DPOSTPROC  -D$(KT) -D$(DT) -D$(IT)
CCFLAGS2 = -O4 -D_FORCE_INLINES -DRODC_ON -DPOSTPROC  -D$(KT) -D$(DT) -D$(IT) 

regex_gpu: main.o topology_function.o general_stdinc.o general_kernels.o
	$(CC)  $(CCFLAGS2) main.o topology_function.o general_stdinc.o general_kernels.o -o regex_gpu

main.o: main.cu
	$(CC) $(CCFLAGS) -c main.cu  

topology_function.o: topology_function.cu
	$(CC) $(CCFLAGS) -c topology_function.cu  

general_stdinc.o: general_stdinc.cpp
	$(CC) $(CCFLAGS) -c general_stdinc.cpp 

general_kernels.o: general_kernels.cu
	$(CC) $(CCFLAGS) -c general_kernels.cu 

clean:
	rm regex_gpu *.o
