#include <stdio.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/matrix.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
R A N K  A L G O R I T H M  R O U T I N E S
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define	MATRIX_FORWARD_ELIMINATION	0
#define	MATRIX_BACKWARD_ELIMINATION	1

int
computeRank(int M, int Q, BitSequence **matrix)
{
	int		i, rank, m=MIN(M,Q);
	
	/* FORWARD APPLICATION OF ELEMENTARY ROW OPERATIONS */ 
	for ( i=0; i<m-1; i++ ) {
		if ( matrix[i][i] == 1 ) 
			perform_elementary_row_operations(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix);
		else { 	/* matrix[i][i] = 0 */
			if ( find_unit_element_and_swap(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix) == 1 ) 
				perform_elementary_row_operations(MATRIX_FORWARD_ELIMINATION, i, M, Q, matrix);
		}
	}

	/* BACKWARD APPLICATION OF ELEMENTARY ROW OPERATIONS */ 
	for ( i=m-1; i>0; i-- ) {
		if ( matrix[i][i] == 1 )
			perform_elementary_row_operations(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix);
		else { 	/* matrix[i][i] = 0 */
			if ( find_unit_element_and_swap(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix) == 1 )
				perform_elementary_row_operations(MATRIX_BACKWARD_ELIMINATION, i, M, Q, matrix);
		}
	} 

	rank = determine_rank(m, M, Q, matrix);

	return rank;
}

void
perform_elementary_row_operations(int flag, int i, int M, int Q, BitSequence **A)
{
	int		j, k;
	
	if ( flag == MATRIX_FORWARD_ELIMINATION ) {
		for ( j=i+1; j<M;  j++ )
			if ( A[j][i] == 1 ) 
				for ( k=i; k<Q; k++ ) 
					A[j][k] = (A[j][k] + A[i][k]) % 2;
	}
	else {
		for ( j=i-1; j>=0;  j-- )
			if ( A[j][i] == 1 )
				for ( k=0; k<Q; k++ )
					A[j][k] = (A[j][k] + A[i][k]) % 2;
	}
}

int
find_unit_element_and_swap(int flag, int i, int M, int Q, BitSequence **A)
{ 
	int		index, row_op=0;
	
	if ( flag == MATRIX_FORWARD_ELIMINATION ) {
		index = i+1;
		while ( (index < M) && (A[index][i] == 0) ) 
			index++;
			if ( index < M )
				row_op = swap_rows(i, index, Q, A);
	}
	else {
		index = i-1;
		while ( (index >= 0) && (A[index][i] == 0) ) 
			index--;
			if ( index >= 0 )
				row_op = swap_rows(i, index, Q, A);
	}
	
	return row_op;
}

int
swap_rows(int i, int index, int Q, BitSequence **A)
{
	int			p;
	BitSequence	temp;
	
	for ( p=0; p<Q; p++ ) {
		temp = A[i][p];
		A[i][p] = A[index][p];
		A[index][p] = temp;
	}
	
	return 1;
}

int
determine_rank(int m, int M, int Q, BitSequence **A)
{
	int		i, j, rank, allZeroes;
	
	/* DETERMINE RANK, THAT IS, COUNT THE NUMBER OF NONZERO ROWS */
	
	rank = m;
	for ( i=0; i<M; i++ ) {
		allZeroes = 1; 
		for ( j=0; j<Q; j++)  {
			if ( A[i][j] == 1 ) {
				allZeroes = 0;
				break;
			}
		}
		if ( allZeroes == 1 )
			rank--;
	} 
	
	return rank;
}

BitSequence**
create_matrix(int M, int Q)
{
	int			i;
	BitSequence	**matrix;
	
	if ( (matrix = (BitSequence **) calloc(M, sizeof(BitSequence *))) == NULL ) {
		printf("ERROR IN FUNCTION create_matrix:  Insufficient memory available.\n");
		
		return NULL;
	}
	else {
		for ( i=0; i<M; i++ ) {
			if ( (matrix[i] = calloc(Q, sizeof(BitSequence))) == NULL ) {
				printf("ERROR IN FUNCTION create_matrix: Insufficient memory for %dx%d matrix.\n", M, M);

				return NULL;
			}
		}
		return matrix;
	}
}

void
def_matrix(int M, int Q, BitSequence **m,int k)
{
	int		i,j;
	
	for ( i=0; i<M; i++ ) 
		for ( j=0; j<Q; j++ )
			m[i][j] = epsilon[k*(M*Q)+j+i*M];
}

void
delete_matrix(int M, BitSequence **matrix)
{
	int		i;

	for ( i=0; i<M; i++ )
		free(matrix[i]);
	free(matrix);
}




void Print_Matrix(unsigned int M[32]){
	int r,c;
	for(r = 0; r < 32; r++)
	{
		for(c = 31; c >= 0; c--)
		{
			if( (M[r]&( 1 << c)) != 0)
				printf("1");
			else printf("0");
		}
		printf("\n");
	}
	printf("\n");
}
void swap(unsigned int* a,unsigned int* b){
	int c;
	c = *a;
	*a = *b;
	*b = c;
}


int Mrank(unsigned int M[32]){
	int i,j,row,rank = 32;
	unsigned int column_value = (1u << 31);
	
	for(row = 0; row < 32; row++)
	{
		i = row;
		while(column_value)
		{
			if( column_value & M[i] )
			{
				//swap(&M[row],&M[i]);
				j = M[row];
				M[row] = M[i];
				M[i] = j;
				break;
			}
			if(i == 31){
				if(--rank < 31) return 30;
				i=row;
				column_value >>= 1;
				
			}
			else ++i;
		}
		for( j = i + 1; j < 32; j++)
		{
			
			/*if( (M[j] & column_value) != 0)	
			{
				M[j] ^= M[row];
			}*/
			//Equivalent to previous if
			M[j] ^= M[row] * ( (M[j] & column_value) / column_value);
		}
		column_value >>= 1;
		//Print_Matrix(M);
	}	
	//Print_Matrix(M);
	
	return rank;
}


