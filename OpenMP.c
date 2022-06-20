#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "OpenMP.h"

Mutant *mutantCreation(Sequence *seq1, Sequence *seq2, double *weights, int typeMutant, int processNum)
{
    int nextOffset = 0, totalOffsets;
    int numOffsets = (seq1->size - seq2->size) + 1;
    if (processNum == 0) // MASTER
         totalOffsets = numOffsets / 2;
    else // SLAVE
    {
	 nextOffset = numOffsets / 2;
	 totalOffsets = numOffsets - nextOffset;
    }
        
    // Allocate memory for ompMutant (return to main)
    Mutant *ompMutant = (Mutant *)malloc(sizeof(Mutant));

    // For compare the first mutant with defult mutant by type (max/min)
    if (typeMutant == 1) //max
        ompMutant->mutantScore = -1111111;
    else //min
        ompMutant->mutantScore = 1111111;

#pragma omp parallel for shared(ompMutant)
    for (int i = 0; i < totalOffsets; i++)  // Run over offsets
    {
        for (int j = 0; j < seq2->size; j++) // Each index of seq2
        {
            Mutant *bestMutantPerThread = createMutant(seq1, seq2, j, i + nextOffset, typeMutant, weights);

#pragma omp critical
            {
                if (ompMutant->mutantScore == -1111111 || ompMutant->mutantScore == 1111111)
                {
                    ompMutant->seq = (char *)malloc(sizeof(char) * seq2->size);
                    strcpy(ompMutant->seq, bestMutantPerThread->seq);
                    ompMutant->size = bestMutantPerThread->size;
                    ompMutant->mutantScore = bestMutantPerThread->mutantScore;
                }
                byMutantType(ompMutant, bestMutantPerThread, typeMutant);
		
            }
            free(bestMutantPerThread->seq);
            free(bestMutantPerThread);
        }
    }
    return ompMutant;
}

Mutant *createMutant(Sequence *seq1, Sequence *seq2, int indexSeq2, int offset, int typeMutant, double *weights)
{
    const char *coservativeGroup[COSER_SIZE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
    const char *semiCoservativeGroup[SEMI_SIZE] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK",
                                                   "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
    const char letters[LETTERS_SIZE] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
                                        'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
                                        'W', 'X', 'Y', 'Z', '-'};
    
    Mutant *myMutant = (Mutant *)malloc(sizeof(Mutant));

    myMutant->seq = NULL;

    for (int i = 0; i < LETTERS_SIZE; i++) // Run over letters array
    {
        char letter = letters[i];

        if (checkPairs(&(seq2->seq[indexSeq2]), letter, coservativeGroup, COSER_SIZE))
        {
            if (myMutant->seq == NULL) // Case of compare the first mutant with defult mutant (seq2)
            {
                myMutant->seq = (char *)malloc(sizeof(char) * seq2->size);

                strcpy(myMutant->seq, seq2->seq);
                myMutant->size = seq2->size;
                myMutant->offset = offset;

                if (typeMutant == 1) //max
                    myMutant->mutantScore = -1111111;
                else //min
                    myMutant->mutantScore = 1111111;
            }

            // Second mutant with new letter to comapre
            Mutant *tempMutant = (Mutant *)malloc(sizeof(Mutant));
            tempMutant->seq = (char *)malloc(sizeof(char) * myMutant->size);
            strncpy(tempMutant->seq, myMutant->seq, myMutant->size);
            tempMutant->size = myMutant->size;
            tempMutant->offset = myMutant->offset;
            tempMutant->seq[indexSeq2] = letter;

            tempMutant->mutantScore = calcMutantScore(seq1, tempMutant, weights, coservativeGroup, COSER_SIZE, semiCoservativeGroup, SEMI_SIZE);
            byMutantType(myMutant, tempMutant, typeMutant);
            free(tempMutant->seq);
            free(tempMutant);
        }
    }
    return myMutant;
}

// Checking if the char from seq2 can be replace according to mutant's rules
int checkPairs(char *charSeq2, char letter, const char **group, int groupSize)
{
    if (*charSeq2 == letter)
        return 0;
    if (checkInGroupForChange(charSeq2, letter, group, groupSize))
        return 1;
    return 0;
}

// 0 -> The char from seq2 and the letter in the same group, 1 -> can be replace
int checkInGroupForChange(char *charSeq2, char letter, const char **group, int groupSize)
{
    for (int i = 0; i < groupSize; i++)
    {
        if (strchr(group[i], *charSeq2) && strchr(group[i], letter))
            return 0;
    }
    return 1;
}

// Calculate Sequence Alignment by using countSignsArr
double calcMutantScore(Sequence *seq1, Mutant *myMutant, double *weights, const char **coservativeGroup, int coserSize, const char **semiCoservativeGroup, int semiSize)
{
    int countSignsArr[] = {0, 0, 0, 0}; // {* , : , . , _}

    for (int i = 0; i < myMutant->size; i++)
    {
        char sign = defineSignsMutate(&(seq1->seq[i]), &(myMutant->seq[i]), myMutant->offset, coservativeGroup, coserSize, semiCoservativeGroup, semiSize);

        if (sign == STAR)
            countSignsArr[0] += 1;
        else if (sign == COLON)
            countSignsArr[1] += 1;
        else if (sign == POINT)
            countSignsArr[2] += 1;
        else //sign == SPACE
            countSignsArr[3] += 1;
    }

    return weights[0] * (double)(countSignsArr[0]) - weights[1] * (double)(countSignsArr[1]) - weights[2] * (double)(countSignsArr[2]) - weights[3] * (double)(countSignsArr[3]);
}

char defineSignsMutate(char *charSeq1, char *mutantChar, int offset, const char **coservativeGroup, int coserSize, const char **semiCoservativeGroup, int semiSize)
{
    if (*(charSeq1 + offset) == *mutantChar)
        return STAR;
    else if (*(charSeq1 + offset) == '-' && *mutantChar != '-')
        return SPACE;
    else if (checkInGroup(charSeq1 + offset, mutantChar, coservativeGroup, coserSize))
        return COLON;
    else if (checkInGroup(charSeq1 + offset, mutantChar, semiCoservativeGroup, semiSize))
        return POINT;
    else
        return SPACE;
}

// Checking if char from seq1 and char from mutant are in a specific group according to pair-wise comparison
int checkInGroup(char *charSeq1, char *mutantChar, const char **group, int groupSize)
{
    for (int i = 0; i < groupSize; i++)
    {
        if (strchr(group[i], *charSeq1) && strchr(group[i], *mutantChar))
            return 1;
    }
    return 0;
}

// Checking which mutant is the maximum or minimum between two mutants
void byMutantType(Mutant *mutant1, Mutant *mutant2, int typeMutant)
{

    if (typeMutant) //max
    {
        if (mutant2->mutantScore > mutant1->mutantScore)
        {
            strncpy(mutant1->seq, mutant2->seq, mutant2->size);
            mutant1->mutantScore = mutant2->mutantScore;
            mutant1->offset = mutant2->offset;
        }
    }
    else // min
    {
        if (mutant2->mutantScore < mutant1->mutantScore)
        {
            strncpy(mutant1->seq, mutant2->seq, mutant2->size);
            mutant1->mutantScore = mutant2->mutantScore;
            mutant1->offset = mutant2->offset;
        }
    }
}
