void setDim(unsigned long int m) {
    unsigned long int k = m / 256;
    int a, b, f;
    unsigned long int i;
    f = sqrt(k);
    i = selectPrime(f);
    a = prime[i + 3];
    b = prime[i];
    dim(a, b);
    printf("2DBF dimensions for aBF: \n%d  %d\n", a, b);
    size = a * b * 64;
}
