#ifndef __OPENMP_H
#define __OPENMP_H

#include "Utils.h"

Mutant *mutantCreation(Sequence *seq1, Sequence *seq2, double *weights, int typeMutant, int processNum);
Mutant *createMutant(Sequence *seq1, Sequence *seq2, int indexSeq2, int offset, int typeMutant, double *weights);
int checkPairs(char *charSeq2, char letter, const char **group, int groupSize);
int checkInGroupForChange(char *charSeq2, char letter, const char **group, int groupSize);
double calcMutantScore(Sequence *seq1, Mutant *myMutant, double *weights, const char **coservativeGroup, int coserSize, const char **semiCoservativeGroup, int semiSize);
char defineSignsMutate(char *charSeq1, char *mutantChar, int offset, const char **coservativeGroup, int coserSize, const char **semiCoservativeGroup, int semiSize);
int checkInGroup(char *charSeq1, char *mutantChar, const char **group, int groupSize);
void byMutantType(Mutant *mutant1, Mutant *mutant2, int typeMutant);


#endif
