#ifndef __UTILS_H
#define __UTILS_H

#define COUNT_SIGNS_SIZE 4
#define COSER_SIZE 9
#define SEMI_SIZE 11
#define LETTERS_SIZE 27
#define SEQ1_SIZE 10000
#define SEQ2_SIZE 5000

extern const char *coservativeGroup[COSER_SIZE];
extern const char *semiCoservativeGroup[SEMI_SIZE];
extern const char letters[LETTERS_SIZE];

// Signs of compareness
#define STAR '*'
#define COLON ':'
#define POINT '.'
#define SPACE ' '

// To decrease the code size
#define CHECK_NULL_INT(value,msg) if(value==NULL){printf(msg);return 0;}
#define CHECK_NULL_VOID(value,msg) if(value==NULL){printf(msg);return;}

typedef struct {
	char *seq; 
	int size;
} Sequence;

typedef struct {
	char *seq; 
	int size;
	int offset;
	double mutantScore; // how much the mutant is similar to the original (seq1)
} Mutant;


#endif 
