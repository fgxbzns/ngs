// filename: foo.cpp

#include <stdio.h>
extern "C"
char* myprint(char *str)
{
    puts(str);
    return str;
}
extern "C"
float add(float a, float b)
{
    return a + b;
}



