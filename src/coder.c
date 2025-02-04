#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "headers/count_digits.h"
#include "headers/buildGF.h"
#include "headers/binaryToDecimal.h"

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
    return 0;  
}
