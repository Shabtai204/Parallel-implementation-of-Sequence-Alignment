#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "FileProcessing.h"
#include "Utils.h"
#include "OpenMP.h"
#include "CudaFuncs.h"

#define WEIGHTS_SIZE 4
#define SCORE_SIZE 7

#define MASTER 0
#define SLAVE 1

int main(int argc, char *argv[])
{
  char *inputFile, *outputFile;
  Sequence *seq1, *seq2;
  double *weights;
  int weightsSize = WEIGHTS_SIZE;
  char *scoreType;
  int typeMutant;
  Mutant *mutantOmp1;
  Mutant *mutantOmp2;
  inputFile = argv[1];
  outputFile = argv[2];

  int my_rank, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Status status;

  if (my_rank == MASTER)
  {
    weights = (double *)malloc(sizeof(double) * WEIGHTS_SIZE);
    CHECK_NULL_INT(weights, "Cannot allocate memory for weights array");

    seq1 = (Sequence *)malloc(sizeof(Sequence));
    CHECK_NULL_INT(seq1, "Cannot allocate memory for seq1");

    seq1->seq = (char *)malloc(sizeof(char) * SEQ1_SIZE);
    CHECK_NULL_INT(seq1->seq, "Cannot allocate memory for seq1");

    seq2 = (Sequence *)malloc(sizeof(Sequence));
    CHECK_NULL_INT(seq2, "Cannot allocate memory for seq2");

    seq2->seq = (char *)malloc(sizeof(char) * SEQ2_SIZE);
    CHECK_NULL_INT(seq2->seq, "Cannot allocate memory for seq2");

    scoreType = (char *)malloc(sizeof(char) * TYPE_SIZE);
    CHECK_NULL_INT(scoreType, "Cannot allocate memory for scoreType");

    mutantOmp1 = (Mutant *)malloc(sizeof(Mutant));
    CHECK_NULL_INT(mutantOmp1, "Cannot allocate memory for mutantOmp1");

    mutantOmp1->seq = (char *)malloc(sizeof(char) * SEQ2_SIZE);
    CHECK_NULL_INT(mutantOmp1->seq, "Cannot allocate memory for  mutantOmp1->seq");

    mutantOmp2 = (Mutant *)malloc(sizeof(Mutant));
    CHECK_NULL_INT(mutantOmp2, "Cannot allocate memory for mutantOmp2");

    mutantOmp2->seq = (char *)malloc(sizeof(char) * SEQ2_SIZE);
    CHECK_NULL_INT(seq1, "Cannot allocate memory for mutantOmp2->seq");

    readFromFile(inputFile, weights, seq1, seq2, scoreType);

    typeMutant = !(strcmp(scoreType, "maximum")); // 1 maximum 0 minimum

    MPI_Send(weights, 4, MPI_DOUBLE, SLAVE, 1, MPI_COMM_WORLD);
    MPI_Send(&seq1->size, 1, MPI_INT, SLAVE, 2, MPI_COMM_WORLD);
    MPI_Send(seq1->seq, seq1->size, MPI_CHAR, SLAVE, 3, MPI_COMM_WORLD);
    MPI_Send(&seq2->size, 1, MPI_INT, SLAVE, 4, MPI_COMM_WORLD);
    MPI_Send(seq2->seq, seq2->size, MPI_CHAR, SLAVE, 5, MPI_COMM_WORLD);
    MPI_Send(&typeMutant, 1, MPI_INT, SLAVE, 6, MPI_COMM_WORLD);

    mutantOmp1 = mutantCreation(seq1, seq2, weights, typeMutant, MASTER);

    MPI_Recv(mutantOmp2->seq, seq2->size, MPI_CHAR, SLAVE, 9, MPI_COMM_WORLD, &status);
    MPI_Recv(&mutantOmp2->offset, 1, MPI_INT, SLAVE, 10, MPI_COMM_WORLD, &status);
    MPI_Recv(&mutantOmp2->mutantScore, 1, MPI_DOUBLE, SLAVE, 11, MPI_COMM_WORLD, &status);

    if (mutantOmp1->mutantScore > mutantOmp2->mutantScore)
      writeToFile(outputFile, mutantOmp1);
    else
       writeToFile(outputFile, mutantOmp2);

    free(mutantOmp2->seq);
    free(mutantOmp2);
    free(mutantOmp1->seq);
    free(mutantOmp1);
    free(scoreType);
    free(seq2->seq);
    free(seq2);
    free(seq1->seq);
    free(seq1);
    free(weights);
  }
  else
  {
    weights = (double *)malloc(sizeof(double *) * WEIGHTS_SIZE);
    CHECK_NULL_INT(seq1, "Cannot allocate memory for weights array");

    seq1 = (Sequence *)malloc(sizeof(Sequence));
    CHECK_NULL_INT(seq1, "Cannot allocate memory for seq1");

    seq1->seq = (char *)malloc(sizeof(char) * SEQ1_SIZE);
    CHECK_NULL_INT(seq1->seq, "Cannot allocate memory for seq1");

    seq2 = (Sequence *)malloc(sizeof(Sequence));
    CHECK_NULL_INT(seq2, "Cannot allocate memory for seq2");

    seq2->seq = (char *)malloc(sizeof(char) * SEQ2_SIZE);
    CHECK_NULL_INT(seq2->seq, "Cannot allocate memory for seq2");

    scoreType = (char *)malloc(sizeof(char) * TYPE_SIZE);
    CHECK_NULL_INT(scoreType, "Cannot allocate memory for scoreType");

    MPI_Recv(weights, 4, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&seq1->size, 1, MPI_INT, MASTER, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(seq1->seq, seq1->size, MPI_CHAR, MASTER, 3, MPI_COMM_WORLD, &status);
    MPI_Recv(&seq2->size, 1, MPI_INT, MASTER, 4, MPI_COMM_WORLD, &status);
    MPI_Recv(seq2->seq, seq2->size, MPI_CHAR, MASTER, 5, MPI_COMM_WORLD, &status);
    MPI_Recv(&typeMutant, 1, MPI_INT, MASTER, 6, MPI_COMM_WORLD, &status);

    mutantOmp2 = mutantCreation(seq1, seq2, weights, typeMutant, SLAVE);

    MPI_Send(mutantOmp2->seq, mutantOmp2->size, MPI_CHAR, MASTER, 9, MPI_COMM_WORLD);
    MPI_Send(&mutantOmp2->offset, 1 , MPI_INT, MASTER, 10, MPI_COMM_WORLD);
    MPI_Send(&mutantOmp2->mutantScore, 1, MPI_DOUBLE, MASTER, 11, MPI_COMM_WORLD);

  }

  MPI_Finalize();
  return 0;
}
