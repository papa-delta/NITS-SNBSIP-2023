
//#include "mask/mask5.h"
//#include "mask/mask6.h"
//#include "mask/mask7.h"
#include "mask/mask8.h"
//#include "mask/mask9.h"
//#include "mask/mask10.h"
//#include "mask/mask12.h"
//#include "mask/mask14.h"
//#include "mask/mask16.h"

//#include "initCBF.h" - commented out initially

static unsigned long int nc=8; //total number of counters
static unsigned long int bc=8; // total number bits per counters

// static int x; - commented because redundant
// static int y;
//int seed[8]={7689571,15485863,104395303,49979693,67867967,141650939,870889,899939};


//Insertion function of KmerCo with only K-mer as input
void _set_(unsigned long int **a,int kmer_len,char *kmer, int k)
{
	unsigned long int h;
	int i,j,pos,loop;
	unsigned long int c;
	for(loop=0;loop<k;loop++){
		h=murmur2(kmer,kmer_len,seed[loop]);
		i=h%x;
		j=h%y;
		pos=h%nc; 
		c=a[i][j]&em[pos];
		c=c>>(bc*(unsigned long int)pos);
		c=c+1UL;		
		if(c==0xFF)
			return;
		c=c<<(bc*pos); 
		a[i][j]=a[i][j]&rm[pos];
		a[i][j]=a[i][j] | c;
	}
}

//if( (c==0x1F && bc==5) || (c==0xFF && bc==8) ||  (c==0x3F && bc==6) || (c==0x7F && bc==7) ||  (c==0x1FF && bc==9)|| (c==0x3FF && bc==10) || (c==0xFFF && bc==12) || (c==0x3FFF && bc==14) ||  (c==0xFFFF && bc==16) ) //Write the condition of current experiment at the beginning for quicker time
//Query function of KmerCo with only K-mer as input
unsigned long int _test_(unsigned long int **a,int kmer_len,char *kmer, int k)
{
	unsigned long int h;
	int i,j,pos,loop;
	unsigned long int c,count[k],min;
	for(loop=0;loop<k;loop++){
		h=murmur2(kmer,kmer_len,seed[loop]);
		i=h%x; 
		j=h%y;
		pos=h%nc; 
		c=a[i][j]&em[pos];
		count[loop]=c>>(bc*(unsigned long int)pos);
		if(count[loop]==0)
			return 0;
	}
	switch(k){
		case 1:
			return count[0];
		case 2:
			if(count[0]<count[1])
				return count[0];
			else
				return count[1];
		default:	
			min=count[0];
			for(loop=1;loop<k;loop++){
				if (min>count[loop])
					min=count[loop];
			}
			return min;
	}
}

//Insertion function of KmerCo with K-mer and reverse complement of K-mer as input
int _set_canonical_(unsigned long int **a,int kmer_len,char *kmer, char *rev_kmer,int k) //h(kmer)<h(rev_kmer)?return 0: return 1
{
	unsigned long int h,h1;
	int i,j,pos,loop,result=0;
	unsigned long int c;
	h=murmur2(kmer,kmer_len,seed[0]); 
   	h1=murmur2(rev_kmer,kmer_len,seed[0]);
	if (h1<h){
		h=h1;
		result=1;
		kmer=rev_kmer;
	}
	for(loop=0;loop<k;loop++){
		i=h%x;
		j=h%y;
		pos=h%nc; 
		c=a[i][j]&em[pos];
		c=c>>(bc*(unsigned long int)pos);
		c=c+1UL;		
		if(c==0xFF)
			return result;
		c=c<<(bc*pos); 
		a[i][j]=a[i][j]&rm[pos];
		a[i][j]=a[i][j] | c;
		if (k>1)
			h=murmur2(kmer,kmer_len,seed[loop+1]); 
	}
	return result;	
}

//if( (c==0x1F && bc==5) || (c==0xFF && bc==8) ||   (c==0x3F && bc==6) || (c==0x7F && bc==7) ||  (c==0x1FF && bc==9)|| (c==0x3FF && bc==10)|| (c==0xFFF && bc==12) || (c==0x3FFF && bc==14) || (c==0xFFFF && bc==16)) //Write the condition of current experiment at the beginning for quicker time
//Query function of KmerCo with K-mer and reverse complement of K-mer as input
unsigned long int _test_canonical_(unsigned long int **a,int kmer_len,char *kmer, char *rev_kmer, int k, int *result)
{
	unsigned long int h,h1;
	int i,j,pos,loop;	
	unsigned long int c,count[k],min;
	*result=0;

	h=murmur2(kmer,kmer_len,seed[0]); 
    	h1=murmur2(rev_kmer,kmer_len,seed[0]);
	if (h1<h){
		h=h1;
		*result=1;
		kmer=rev_kmer;
	}		
	for(loop=0;loop<k;loop++){
		i=h%x; 
		j=h%y;
		pos=h%nc; 
		c=a[i][j]&em[pos];
		count[loop]=c>>(bc*(unsigned long int)pos);
		if(count[loop]==0)
			return 0;
		if (k>1)
			h=murmur2(kmer,kmer_len,seed[loop+1]); 
	}
	switch(k){
		case 1:
			return count[0];
		case 2:
			if(count[0]<count[1])
				return count[0];
			else
				return count[1];
		default:	
			min=count[0];
			for(loop=1;loop<k;loop++){
				if (min>count[loop])
					min=count[loop];
			}
			return min;
	}
}
