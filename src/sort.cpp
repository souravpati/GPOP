#include "../include/sort.h"

unsigned int computeMean (unsigned int a, unsigned int b)
{
    unsigned int meanVal = (a >> 1) + (b >> 1) + ((a & LSB_MASK) + (b & LSB_MASK))/2;
    return meanVal;
}

