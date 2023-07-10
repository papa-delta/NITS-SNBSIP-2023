
void setDim(unsigned long int m)
{
	unsigned long int k=m/256; 
	int a,b,c,d,e,f;
	unsigned long int i;
	f=sqrt(k);
	i=selectPrime(f);
	a=prime[i+3];
	b=prime[i];
	//c=prime[i/2-3];
	//d=prime[i/2+3];
	dim(a,b);
	//dim2(c,d);
	printf("2DBF dimensions for aBF: \n%d  %d\n",a,b);
	//printf("2DBF dimensions for bBF: \n%d  %d\n",c,d);
	size=a*b*64;

}
