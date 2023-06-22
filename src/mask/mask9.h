//9-bit counter mask

unsigned long int rm[7]={
    0xfffffffffffffE00,
    0xfffffffffffC01ff,
    0xfffffffff803ffff,
    0xfffffff007ffffff,
    0xffffE00fffffffff,
    0xffC01fffffffffff,
    0x803fffffffffffff};
unsigned long int em[7]={
    0x00000000000001ff,
    0x000000000003fE00,
    0x0000000007fC0000,
    0x0000000ff8000000,
    0x00001ff000000000,
    0x003fE00000000000,
    0x7fC0000000000000};
