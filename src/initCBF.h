
unsigned long int size=0;

unsigned long int selectPrime(unsigned long int k)
{
	unsigned long int i;
	for(i=1;i<total_prime;i++)
	{
		if(prime[i]>k)
			return i;
	}
}

double error(unsigned long int m, unsigned long int n)
{
	return pow((1-exp(-2*n/m)),2);
}
unsigned long int memory(unsigned long int n, double err)
{
	return (unsigned long int)(-(n*log(err))/pow(log(2),2));
}
unsigned long int number(unsigned long int m,double err)
{
	return (unsigned long int)(-(m*pow(log(2),2))/log(err));
}

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
