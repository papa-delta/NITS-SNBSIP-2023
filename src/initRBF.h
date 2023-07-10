
void setDim_rBF(unsigned long int m)
{
	unsigned long int k=m/(2*64); 
	int a,b,c,d,e,f;
	unsigned long int i;

	f=sqrt(k);
	i=selectPrime(f);
	//int j=(i/1.5);
	a=prime[i/2+3];
	b=prime[i/2-3];
	//c=prime[i-3];
	//d=prime[i+3];
	//e=prime[i/3+3];
	//f=prime[i/3-3];
	dim(a,b);
	//dim1(c,d);
	//dim2(e,f);
	 
	printf("2DBF dimensions: \n%d  %d\n",a,b);
	size+=(a*b)*64;

}
