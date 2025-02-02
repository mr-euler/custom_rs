#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int count_digits(int i) { //считаем количество цифр
    struct {
            int max;
            int count;
    } static digits[32] = {
            { 9, 1 }, { 9, 1 }, { 9, 1 }, { 9, 1 },
            { 99, 2 }, { 99, 2 }, { 99, 2 },
            { 999, 3 }, { 999, 3 }, { 999, 3 },
            { 9999, 4 }, { 9999, 4 }, { 9999, 4 }, { 9999, 4 },
            { 99999, 5 }, { 99999, 5 }, { 99999, 5 },
            { 999999, 6 }, { 999999, 6 }, { 999999, 6 },
            { 9999999, 7 }, { 9999999, 7 }, { 9999999, 7 }, { 9999999, 7 },
            { 99999999, 8 }, { 99999999, 8 }, { 99999999, 8 },
            { 999999999, 9 }, { 999999999, 9 }, { 999999999, 9 },
            { 2147483647, 10 }, { 2147483647, 10 }
    };
        register const int z = 0;
        register unsigned log2;
        if (i < 0) i = -i;
        __asm__ __volatile__ (
                "bsr %1, %0;"  \
                "cmovz %2, %0;"\
                : "=r" (log2)  \
                : "rm" (i), "r"(z));
        return digits[log2].count + ( i > digits[log2].max );
}

void buildGF(int *gf, short int n, int polinom)//
{
    for ( int i=1;i<n;i++)
    {
        gf[i]=gf[i-1]>>1;
        if((gf[i-1]%2)==1)
            {
                gf[i]^=polinom;
            }
    }
}

int binaryToDecimal(int num,  int *p) 
{ 
    int count=0;
    int dec_value = 0; 
    int base = 1;
    int last_digit = num % 10; 
    while (num) { 
      //3   int last_digit = num % 10; 
        num /= 10; 
        dec_value += last_digit * base; 
        base *= 2;
        last_digit = num % 10;
        count++;

    } 
    
    *p=count;
    return dec_value; 
}

int main(int args_number, char **args_value) // вводить полином, начиная с младшей стемени. x4+x+1 => ./codec -p 1101
{
    int poly;
    int poly_lengh=5; // контролирует размер полинома
    int n, k, p, t;

    printf("Enter polinom. Example: x+1+x4 => 11001\n");
    scanf("%d",&poly);
    
    if (count_digits(poly)!=poly_lengh)
        {
            printf("polynom should be %d digits long \n", poly_lengh);
            return 0;
        }
    
    printf("Entern info bits lengh \"k\"\n");    
    scanf("%d",&k); //101110001

    poly=binaryToDecimal(poly/10, &p); // получаем вычет в 10й форме

    printf("poly: %d \n", poly);

    n=1<<p;
    n--; // размер поля

    int gf[n];
    int gf_size=sizeof(gf)/sizeof(gf[0]);
    
    gf[0]=1<<p-1;

    buildGF(&gf[0],gf_size,poly);

    /*for(int i=0;i<15;i++)
    {
        printf("gf[%d] = %d \n",i, gf[i]);
    }
    */
    t=(n-k)/2; // количество гарантированно исправляемых ошибок
 
    printf("k: %d \n", t);
    printf("p: %d \n", p);

    return 0;  
}
