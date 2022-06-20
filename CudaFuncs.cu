/*
#include <cuda_runtime.h>
#include <string.h>
#include <stdio.h>

Mutant *setGPUforMutantCreation(Sequence *seq1, Sequence *seq2, int gpuSize, double *weights, int typeMutant)
{
    char *d_seq1, *d_seq2;
    double *d_weights;
    Mutant *gpuMutant;
    Mutant *tempMutant;

    int numOfBlocks = gpuSize;                      // each block 1 offset
    int NumOfThreads = (gpuSize * seq2->size) * 32; // change letters per index

    // Set the allocations for the device
    cudaMalloc((void **)&d_seq1, (seq1->size) * sizeof(char));
    cudaMalloc((void **)&d_seq2, (seq2->size) * sizeof(char));
    cudaMalloc((void **)&d_weights, 4 * sizeof(double));
    cudaMalloc((void **)&gpuMutant, (seq2->size) * sizeof(Mutant));
    cudaMalloc((void **)&tempMutant, (seq2->size) * sizeof(Mutant));
    cudaMemcpy(d_seq1, seq1->seq, (seq1->size) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, seq2->seq, (seq2->size) * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_weights, weights, 4 * sizeof(double), cudaMemcpyHostToDevice);
    // Launch Kernel
    mutantCreation<<<numOfBlocks, NumOfThreads>>>(d_seq1, d_seq2, size, gpuSize, d_weights, typeMutant, tempMutant);
    // Get the result from the GPU.
    cudaMemcpy(gpuMutant, tempMutant, (tempMutant->size) * sizeof(Mutant), cudaMemcpyDeviceToHost);
    cudaFree(tempMutant->seq);
    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(d_weights);
    return gpuMutant;
}
// 1 thread
__global__ void *mutantCreation(Sequence *seq1, Sequence *d_seq2, int gpuSize, double *weights, int typeMutant, Mutant *tempMutant)
{
    int numOfThreadsPerBlock = seq2->size;
    int myIndex = threadIdx.x + blockIdx.x * numOfThreadsPerBlock;

    __shared__ Mutant *tempMutant;
    cudaMalloc((void **)&tempMutant, (seq2->size) * sizeof(Mutant));

    // To compare the first mutant with defult mutant by type (max/min)
    if (typeMutant == 1) //max
        tempMutant->mutantScore = -1111111;
    else //min
        tempMutant->mutantScore = 1111111;

    Mutant *bestMutantPerThread = createMutant(seq1, seq2,  threadIdx.x,  blockIdx.x, typeMutant, weights);

    if (tempMutant->mutantScore == -1111111 || tempMutant->mutantScore == 1111111)
    {
        cudaMalloc((void **)&(tempMutant->seq), (seq2->size) * sizeof(char)); //tempMutant->seq = (char *)malloc(sizeof(char) * seq2->size);
        if (tempMutant->seq == NULL)
        {
            printf("Cannot allocate memory for myMutant");
        }
        cudaMemcpy(tempMutant->seq, bestMutantPerThread->seq, (bestMutantPerThread->size) * sizeof(Mutant), cudaMemcpyDeviceToHost);
        tempMutant->size = bestMutantPerThread->size;
        tempMutant->mutantScore = bestMutantPerThread->mutantScore;
    }
    tempMutant->offset = blockIdx.x + threadIdx.x;
    byMutantType(tempMutant, bestMutantPerThread, typeMutant);
    cudaFree(bestMutantPerThread->seq);
    cudaFree(bestMutantPerThread->seq);

    __syncthreads(); // wait for all threads to finish
}

__device__ Mutant *createMutant(Sequence *seq1, Sequence *seq2, int indexSeq2, int offset, int typeMutant, double *weights)
{
    const char *coservativeGroup[COSER_SIZE] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
    const char *semiCoservativeGroup[SEMI_SIZE] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK",
                                                   "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
    const char letters[LETTERS_SIZE] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
                                        'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
                                        'W', 'X', 'Y', 'Z', '-'};
    Mutant *myMutant;
    cudaMalloc((void **)&myMutant, sizeof(Mutant)); //Mutant *myMutant = (Mutant *)malloc(sizeof(Mutant));
    if (myMutant == NULL)
    {
        printf("Cannot allocate memory for myMutant");
        return NULL;
    }
    myMutant->seq = NULL;

    for (int i = 0; i < LETTERS_SIZE; i++)
    {
        char letter = letters[i];

        if (checkPairs(&(seq2->seq[indexSeq2]), letter, coservativeGroup, COSER_SIZE))
        {
            if (myMutant->seq == NULL) // Case of compare the first mutant with defult mutant (seq2)
            {
                cudaMalloc((void **)&myMutant->seq, (seq2->size) * sizeof(Mutant)); //myMutant->seq = (char *)malloc(sizeof(char) * seq2->size);

                if (myMutant->seq == NULL)
                {
                    printf("Cannot allocate memory for myMutant");
                    return NULL;
                }
                cudaMemcpy(myMutant->seq, seq2->seq, (seq2->size) * sizeof(char), cudaMemcpyDeviceToHost);
                myMutant->size = seq2->size;
                myMutant->offset = offset;

                if (typeMutant == 1) //max
                    myMutant->mutantScore = -1111111;
                else //min
                    myMutant->mutantScore = 1111111;
            }

            // Second mutant to comapre
            cudaMalloc((void **)&tempMutant, sizeof(Mutant)); //Mutant *tempMutant = (Mutant *)malloc(sizeof(Mutant));

            if (tempMutant == NULL)
            {
                printf("Cannot allocate memory for tempMutant");
                return NULL;
            }

            cudaMalloc((void **)&tempMutant->seq, (myMutant->size) * sizeof(char)); //tempMutant->seq = (char *)malloc(sizeof(char) * myMutant->size);
            if (myMutant->seq == NULL)
            {
                printf("Cannot allocate memory for tempMutant");
                return NULL;
            }
            cudaMemcpy(tempMutant->seq, myMutant->seq, (myMutant->size) * sizeof(char), cudaMemcpyDeviceToHost);
            tempMutant->size = myMutant->size;
            tempMutant->offset = myMutant->offset;
            tempMutant->seq[indexSeq2] = letter;
            tempMutant->mutantScore = calcMutantScore(seq1, tempMutant, weights, coservativeGroup, COSER_SIZE, semiCoservativeGroup, SEMI_SIZE);
            byMutantType(myMutant, tempMutant, typeMutant);
            cudaFree(tempMutant->seq);
            cudaFree(tempMutant);
        }
    }
    return myMutant;
}

__device__ int checkPairs(char *charSeq2, char letter, const char **group, int groupSize)
{
    if (*charSeq2 == letter)
        return 0;
    if (checkInGroupForChange(charSeq2, letter, group, groupSize))
        return 1;
    return 0;
}

__device__ int checkInGroupForChange(char *charSeq2, char letter, const char **group, int groupSize)
{
    for (int i = 0; i < groupSize; i++)
    {
        if (_strchr(group[i], *charSeq2) && _strchr(group[i], letter))
            return 0;
    }
    return 1;
}

__device__ double calcMutantScore(Sequence *seq1, Mutant *myMutant, double *weights, const char **coservativeGroup, int coserSize, const char **semiCoservativeGroup, int semiSize)
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

__device__ char defineSignsMutate(char *charSeq1, char *mutantChar, int offset, const char **coservativeGroup, int coserSize, const char **semiCoservativeGroup, int semiSize)
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

__device__ int checkInGroup(char *charSeq1, char *mutantChar, const char **group, int groupSize)
{
    for (int i = 0; i < groupSize; i++)
    {
        if (_strchr(group[i], *charSeq1) != NULL && _strchr(group[i], *mutantChar))
            return 1;
    }
    return 0;
}

__device__ void byMutantType(Mutant *mutant1, Mutant *mutant2, int typeMutant)
{
    if (typeMutant) //max
    {
        if (mutant2->mutantScore > mutant1->mutantScore)
        {
            cudaMemcpy(mutant1->seq, mutant2->seq, (mutant2->size) * sizeof(char), cudaMemcpyDeviceToHost);
            mutant1->mutantScore = mutant2->mutantScore;
            mutant1->offset = mutant2->offset;
        }
    }
    else // min
    {
        if (mutant2->mutantScore < mutant1->mutantScore)
        {
            cudaMemcpy(mutant1->seq, mutant2->seq, (mutant2->size) * sizeof(char), cudaMemcpyDeviceToHost);
            mutant1->mutantScore = mutant2->mutantScore;
            mutant1->offset = mutant2->offset;
        }
    }
}

// Checks if char c is in string s on GPU device
__device__ char *_strchr(const char *s, int c) {
	while (*s != (char) c)
		if (!*s++)
			return 0;
	return (char *) s;
}

*/
