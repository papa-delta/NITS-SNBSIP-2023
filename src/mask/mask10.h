//10-bit counter mask

unsigned long int rm[6]={
    0xfffffffffffffC00,
    0xfffffffffff003ff,
    0xffffffffC00fffff,
    0xffffff003fffffff,
    0xfffC00ffffffffff,
    0xf003ffffffffffff};
unsigned long int em[6]={
    0x00000000000003ff,
    0x00000000000ffC00,
    0x000000003ff00000,
    0x000000ffC0000000,
    0x0003ff0000000000,
    0x0ffC000000000000};
