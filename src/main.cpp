#include <stdio.h>
#include <RLC.h>
#include <matplotlib.hpp>
void swap(int & a, int & b)
{
    int temp;
    temp = b;
    b = a;
    a = temp;
}

int main(int argc, char const *argv[])
{
    RLCBranch rlc;
    
    int a; 
    int b;
    a = 1; //0x00-0x03
    b = 10;
    rlc.R=1;
    rlc.L=3;
    rlc.C=10;
    RLCBranch * a_pointer = &rlc;
    Real * temp = (Real *)(a_pointer);
    (*temp)++;
    for(int i=0;i<10;i++){
        a++;
    }
    
    printf("%f %d\n",temp[1] , b);
    swap(a, b);
    printf("%d %d\n", a, b);
    printf("Hello World!\n");
    return 0;
}
