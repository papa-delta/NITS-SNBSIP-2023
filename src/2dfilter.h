//#include<stdio.h>
//#include<stdlib.h>
//#include"murmur.h"

//int seed=98899;
//int seed1=10007; 
//int x=103,y=107,z=109; //diemnsion of the 2D filter.
//int x=9973,y=10007;//,z=137;  //dimension of the 2D filter. A prime number is a candidate for good hashing.
//int x=5647,y=6151;
//int x2=1129,y2=1123;
//int x2=179, y2=197;//275KB

//static int x,y; - commented because redundant


//#include "common.h"

static int bits=61;

void _set_(unsigned long int **a,unsigned long int h)
{
	int i,j,pos;
	i=h%x;
	j=h%y;
	pos=h%bits; 
	a[i][j]=a[i][j]|(1UL<<pos);
}


int _test_(unsigned long int **a,unsigned long int h)
{
	int i,j,pos;
	int flag;
	i=h%x; 
	j=h%y;
	pos=h%bits;
	return ((a[i][j]&(1UL<<pos))>>pos);
}

void _del_(unsigned long int **a,unsigned long int h)
{
	int i,j,pos;
	unsigned long int p;
	int flag;
	i=h%x; 
	j=h%y;
	pos=h%bits;
	p=(1UL<<pos);
	if(p==(a[i][j]&(1UL<<pos)))
		a[i][j]=a[i][j]^p;
}



