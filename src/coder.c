

#include <stdio.h>
#include "gflib/gflib.h"
#include "gflib/polinom.h"
#include "gflib/gen_polinom.h"

#define gf_get gf_get_by_degree

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

//void poly_division()


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
    //form_polinom=10011;
    form_polinom=1011;
    form_polinom = decimal_to_binary(form_polinom); // Преобразование 1011 => 11
    printf("\n"); // Отступ

    // Количество информационных символов
    int k;
    printf("Введите количество информационных символов.\n");
    printf(">>");
    //scanf("%d", &k);
    k=3;
    //k=9;
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


   // int data[]={12,2,8,1,4,5,10,10,4}; // мой вариант
    //int data[]={10,13,5,8,9,2,7,7,7}; // пример от Гали
 //  int data[]={13,9,7,14,5,3,2,8,14}; // хуйня со студфайлов
    int data[]={4,-1,3}; // полином из методички Владимирова
    polinom_t *info_poly=polinom_init(gf,k);
    for(int i=0;i<k;i++)
    {
       polinom_set(info_poly, i, gf_get_by_id(gf,data[i]+1));
       //polinom_set(info_poly, i, gf_get_by_degree(gf,data[i]));
    }

    polinom_t *encoded=polinom_init(gf,n);
    for(int i=r;i<n;i++)
    {
       polinom_set(encoded, i, info_poly->data[i-r]);
    }

    int temp_divisor[n]; // массив для хранения остатков от деления
    for (int i=0;i<n;i++) temp_divisor[i]=encoded->data[i];

    int multiplier=temp_divisor[n-1];
    int division_result[gen_polinom->degree];

    for(int i=0;i<=n-gen_polinom->degree;i++) // до тех пор пока степень делимого не будет равна степени генераторного полинова, включая её
    {  
        for(int j=0;j<gen_polinom->degree-1;j++)
        {
            // умножаем ген. полином на элемент при старшей степени и складываем его с делимым/остатком
            division_result[j]=gf_add(gf_mult(gf,gen_polinom->data[j],multiplier),temp_divisor[j+k-1-i]);
            temp_divisor[n-gen_polinom->degree-i+j]=division_result[j];
        }
        multiplier=division_result[gen_polinom->degree-2];
    }

    for(int i=0;i<r;i++) encoded->data[i]=temp_divisor[i];
    for(int i=0;i < n;i++) printf("encoded%d = %d\n",i,gf->rev_table[encoded->data[i]]-1);
    printf("\n"); // Отступ

    //for (int i=0;i<n;i++) temp_divisor[i]=encoded->data[i];
    //temp_divisor[0]=gf_get(gf,2); // Раскоментировать, чтобы получился полином Владимирова. Результат проверен, сходится с листочком.
    //temp_divisor[1]=gf_get(gf,2);
    //temp_divisor[2]=gf_get(gf,1);
    //temp_divisor[3]=gf_get(gf,4);
    
  // Освобождение памяти, выделенной под данные
    polinom_free(gen_polinom);
    gf_free(gf);

    return 0;
}