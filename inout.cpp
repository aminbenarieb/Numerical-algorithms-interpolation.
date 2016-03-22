#include "inout.h"
#include <stdio.h>

void input_limits(double &A, double &B, int &parts)
{
	printf("Введите границы таблицы [ a < b ]: ");
	do
		scanf("%lf %lf", &A, &B);
	while (!(A<B));

	printf("Введите число строк в таблицы: ");
	do
		scanf("%d", &parts);
    
	while (!(parts>0));
}

void print_table(double **table, int rows, int cols)
{
	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<cols; j++)
		{
			if (j==0)
				printf("%3.5lf", table[j][i]);
			else
				printf(" | %.5lf", table[j][i]);
		}
		printf("\n");
	}
}


void print_row(double *row, int n)
{
	for (int i=0; i<n; i++)
	{
		if (i!=0)
			printf(" | ");
		printf("%3.2lf", row[i]);
	}
	printf("\n");
}
void print_col(double *col, int n)
{
	for (int i=0; i<n; i++)
	{
		printf("%3.2lf\n", col[i]);
	}
	printf("\n");
}