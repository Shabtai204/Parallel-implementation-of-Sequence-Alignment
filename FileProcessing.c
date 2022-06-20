#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "FileProcessing.h"

void readFromFile(char *fileName, double *weights, Sequence *seq1, Sequence *seq2, char *scoreType)
{
    FILE *f = fopen(fileName, "r");

    if (f == NULL)
    {
        fprintf(stderr, "Error, file %s does not exist...\n", fileName);
        return;
    }

    fscanf(f, "%lf %lf %lf %lf", &weights[0], &weights[1], &weights[2], &weights[3]);
    CHECK_NULL_VOID(weights, "No weights in the file.");
    fscanf(f, "%s\n", seq1->seq);
    CHECK_NULL_VOID(seq1, "No seq1 in the file.");
    fscanf(f, "%s\n", seq2->seq);
    CHECK_NULL_VOID(seq2, "No seq2 in the file.");
    fscanf(f, "%s\n", scoreType);
    CHECK_NULL_VOID(scoreType, "No scoreType in the file.");
    
    seq1->size = strlen(seq1->seq);
    seq2->size = strlen(seq2->seq);
    fclose(f);
}

void writeToFile(char *fileName, Mutant *myMutant)
{
    FILE *f = fopen(fileName, "w");

    if (f == NULL)
    {
        fprintf(stderr, "Error, file %s does not exist...\n", fileName);
        return;
    }

    fprintf(f, "%s\n", myMutant->seq);
    fprintf(f, "Offset: %d, Alignment Score: %.2f", myMutant->offset, myMutant->mutantScore);

    fclose(f);
}
