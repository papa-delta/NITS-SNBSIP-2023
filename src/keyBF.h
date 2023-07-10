

void insert_keyBF(unsigned long int **b2, char *buff, int i)
{
	
	unsigned long int h1,h2;
	//,h3,h4,h5,h6,h7;

	h1=murmur2(buff,i,seed[0]);
	h2=murmur2(buff,i,seed[1]);
	//h3=murmur2(buff,i,98899); //should we change 98899 to 104395303 (seed[2])?
	//h4=murmur2(buff,i,seed4);
	//h5=murmur2(buff,i,seed5);
	
	_setrBF_(b2,h1);
	_setrBF_(b2,h2);
	//_setrBF_(b2,h3);
	//_set_(b2,h4);
	//_set_(b2,h5);
	
}


int lookup_keyBF(unsigned long int **b2, char *buff, int i)
{
	unsigned long int h1,h2;
	
	//,h3,h4,h5,h6,h7;

	h1=murmur2(buff,i,seed[0]);
	h2=murmur2(buff,i,seed[1]);
	//h3=murmur2(buff,i,98899); //should we change 98899 to 104395303 (seed[2])?
	//h4=murmur2(buff,i,seed4);
	//h5=murmur2(buff,i,seed5);	
	if(_testrBF_(b2,h1)==1)
		if(_testrBF_(b2,h2)==1)
			//if(_testrBF_(b2,h3)==1)
				//if(_test_(b2,h4)==1)
					//if(_test_(b2,h5)==1)
						return 1;
	return 0;									
}

void free_keyBF(unsigned long int **b2)
{
	_free_(b2);
}

