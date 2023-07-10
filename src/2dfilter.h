
static int bits=61;

void _setrBF_(unsigned long int **a,unsigned long int h)
{
	int i,j,pos;
	i=h%x;
	j=h%y;
	pos=h%bits; 
	a[i][j]=a[i][j]|(1UL<<pos);
}


int _testrBF_(unsigned long int **a,unsigned long int h)
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



