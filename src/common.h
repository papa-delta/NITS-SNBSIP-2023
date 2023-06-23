//-------might be redundant, might remove later
#include <stdio.h> 
#include <stdlib.h>
// ----------------------------------------------


static int x;
static int y;

int seed[8]={7689571,15485863,104395303,49979693,67867967,141650939,870889,899939};


void dim(int p, int q)
{
	x=p;
	y=q;
}


unsigned long int **allocate()
{
	int i,j,k;
	unsigned long int **a=(unsigned long int **)malloc(x*sizeof(unsigned long int *));
	if(a==NULL)
	{
		printf("Unable to allocate!\n");
		return NULL;
	}
	
	for(i=0;i<x;i++)
	{
		a[i]=(unsigned long int *)malloc(y*sizeof(unsigned long int));
		if(a[i]==NULL)
		{
			printf("Unable to allocate!\n");
			return NULL;
		}
	}
	for(i=0;i<x;i++)
		for(j=0;j<y;j++)
			a[i][j]=0;
	printf("\nAllocated and Initilized 2DBF Successfully...\n");
	return a;
}


void _free_(unsigned long int **a)
{
	free(a);
	printf("\nMemory freed successfully...\n");
}
