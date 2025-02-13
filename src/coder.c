

#include <stdio.h>
#include "gflib/gflib.h"
#include "gflib/polinom.h"
#include "gflib/gen_polinom.h"

int count_bits(int num) {
    int count = 0;
    while (num) {
        num = num >> 1;
        count++;
    }
    return count;
}

int decimal_to_binary(int num) {
    int decimal = 0;
    int mask = 1;
    while (num) {
        decimal += mask * (num % 10);
        num /= 10;
        mask = mask << 1;
    }
    return decimal;
}

int main() {

    // Параметры кода РС

    // Константы

    int b = 0; // Параметр кодирования

    // Динамические параметры

    // Полином
    // int form_polinom = 0b1011; // x^3 + x + 1
    int form_polinom;
    printf("Введите образующий полином в 10-тичной форме.\n");
    printf("Пример: x^3+x+1 => 1011\n");
    printf(">>");
    //scanf("%d", &form_polinom);
    form_polinom=1011;
    form_polinom = decimal_to_binary(form_polinom); // Преобразование 1011 => 11
    printf("\n"); // Отступ

    // Количество информационных символов
    int k;
    printf("Введите количество информационных символов.\n");
    printf(">>");
    //scanf("%d", &k);
    k=3;
    printf("\n"); // Отступ


    // Остальные параметры

    int m = count_bits(form_polinom)-1; // GF(q^m)

    int n = (1 << m)-1; // Количество кодовых символов
    // Подрязумевается, что n+1 == 2^m

    int r = n - k; // Количество исправляющих символов
    int t = r / 2; // Количество возможных для исправления символов 

    // Формирование поля Галуа
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("gf build error: polinom is not primitive\n");
        return 1;
    }
    // gf_print(gf); // Отобразим элементы поля Галуа

    // Формирование порождающего полинома
    polinom_t *gen_polinom = generating_polinom(gf, b, t);
    printf("Порождающий полином:\n");
    polinom_print(gen_polinom);


    int data[]={4,0,3};

    polinom_t *info_poly=polinom_init(k,gf);
    
    for(int i=0;i<k;i++)
    {
       polinom_set(info_poly, i, gf_get(gf,data[i]+1));
    }

    polinom_t *encoded=polinom_init(n,gf);
    //мы здесь

    int coded[n];
    int temp_array[n];
    for (int i=0;i<n;i++) coded[i]=0;
    for (int i=0;i<n;i++) temp_array[i]=0;

    for (int i=n-1;i>=n-k;i--) coded[i]=data[(n-i-1)];


  // for(int i=0;i < n;i++) printf("data%d = %d\n",i,coded[i]);

   printf("\n"); // Отступ

   for(int i=0;i < gen_polinom->degree;i++)
   {
       printf("gen_poly%d = %d\n",i,gf->rev_table[gen_polinom->data[i]]-1);
   }

    int index=n-1;

    int value;
    int temp[n];
    for(int i=0;i < n;i++) temp[i]=coded[i];
    for(int i=0;i < n;i++) printf("temp%d = %d\n",i,temp[i]);
    printf("\n"); // Отступ
    while(index>=gen_polinom->degree)
    {
        value=temp[index];
        printf("value = %d\n",value);
        int temp_devision[index-1]; 
        for(int i=0;i<index-1;i++)
        {
            temp_devision[i]=gf_add(gf_mult(gen_polinom->data[i], gf_get(gf,value+1),gf),temp[i]);
            //temp_devision[i]=gf->rev_table[gf->table[(gf->rev_table[gen_polinom->data[i]]-1+value-2) % (gf->total_quantity-1)+1]^gf->table[temp[i]]]-1;
            //printf("i=%d, poly=%d\n",i,gf->rev_table[gen_polinom->data[i]]-1+value-2);
            printf("res_%d = %d\n",i,gf->rev_table[temp_devision[i]]-1);
        }
        index=0;
    }
    printf("\n"); // Отступ



    // Освобождение памяти, выделенной под данные
    polinom_free(gen_polinom);
    gf_free(gf);

    return 0;
}