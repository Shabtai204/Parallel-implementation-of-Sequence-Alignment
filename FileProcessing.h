#ifndef __FILEPROCESSING_H
#define __FILEPROCESSING_H

#include "Utils.h"

#define SEQ1_SIZE 10000
#define SEQ2_SIZE 5000
#define TYPE_SIZE 7

void readFromFile(char *fileName, double *weights, Sequence *seq1, Sequence *seq2, char *scoreType);
void writeToFile(char *fileName, Mutant *myMutant);

#endif
