#include "keccak_permutation.h"
uint64_t vector_to_bin64(uint8_t* x)
{
    uint64_t result=0;
    for(int i=63;i>=0;i--)
    {
        result+=x[i]*(uint64_t)(pow(2,63-i));
    }
    return result;
}
int main()
{
    FILE *f =fopen("result.txt","r");
    char a[2];
    char b[13];
    int c=fscanf(f,"%s%s",a,b);
    printf("%s,%s\n",a,b);
    if(b[0]='S')
    {
        printf("111");
    }
}