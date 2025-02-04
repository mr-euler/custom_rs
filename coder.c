#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "count_digits.h"
#include "buildGF.h"
#include "binaryToDecimal.h"

int main(int args_number, char **args_value) // вводить полином, начиная с младшей стемени. x4+x+1 => ./codec -p 1101
{
    int poly;
    //int poly_lengh=5; // контролирует размер полинома
    int n, k, p, t;

    //printf("Enter polinom. Example: x+1+x4 => 11001\n");
    //scanf("%d",&poly);
    poly=11001;
    
    //printf("Entern info bits lengh \"k\"\n");    
    //scanf("%d",&k); //101110001

    k=9;
    poly=binaryToDecimal(poly/10, &p); // получаем вычет в 10й форме p принимает значение длины вычета

    printf("poly: %d \n", poly);

    n=1<<p;
    n--; // размер поля

    int gf[n];
    int reverse_gf[n+1]; // значения поля в степени // поле на 1 больше в размере, т.к. первый элемент -1 (степень числа 0)
	int gf_size=n;
	int reverse_gf_size=n+1;
    
    gf[0]=1<<p-1;

    buildGF(&gf[0],&reverse_gf[0],gf_size,poly); // сразу строим поле галуа и обратное поле 

    int generator[n-k+1];
    for(int i=0;i<n-k+1;i++)
    {
        generator[i]=0;
    }

	for(int i=0;i<n;i++)
	{
		printf("GF[%d] = %d\n",i,gf[i]);
	}
	printf("\n");
	for(int i=0;i<n+1;i++)
	{
		printf("reverse_GF[%d] = %d\n",i,reverse_gf[i]);
	}

    return 0;  
}
