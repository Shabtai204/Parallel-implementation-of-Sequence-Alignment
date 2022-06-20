build:
	mpicxx -fopenmp -c main.c -o main.o
	mpicxx -fopenmp -c FileProcessing.c -o FileProcessing.o
	mpicxx -fopenmp -c OpenMP.c -o OpenMP.o
	mpicxx -fopenmp -o runme main.o FileProcessing.o OpenMP.o -ldl -lrt

clean:
	rm -f *.o ./runme

run:
	mpiexec -np 2 ./runme input.txt output.txt

runOn2:
	mpiexec -np 2 -machinefile  mf  -map-by  node  ./runme input.txt output.txt

