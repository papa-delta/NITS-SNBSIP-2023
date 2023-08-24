void setDim_rBF(unsigned long int m) {
    unsigned long int k = m / (2 * 64);
    int a, b, f;
    unsigned long int i;

    f = sqrt(k);
    i = selectPrime(f);
    a = prime[i / 2 + 3];
    b = prime[i / 2 - 3];
    dim(a, b);

    printf("2DBF dimensions: \n%d  %d\n", a, b);
    size += (a * b) * 64;
}
