
//8-bit counter mask

unsigned long int rm[8]={
    0xffffffffffffff00,
    0xffffffffffff00ff,
    0xffffffffff00ffff,
    0xffffffff00ffffff,
    0xffffff00ffffffff,
    0xffff00ffffffffff,
    0xff00ffffffffffff,
    0x00ffffffffffffff};
unsigned long int em[8]={
    0x00000000000000ff,
    0x000000000000ff00,
    0x0000000000ff0000,
    0x00000000ff000000,
    0x000000ff00000000,
    0x0000ff0000000000,
    0x00ff000000000000,
    0xff00000000000000};

