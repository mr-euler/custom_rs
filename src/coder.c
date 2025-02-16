

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

    int b = 1; // Параметр кодирования

    // Динамические параметры

    // Полином
    // int form_polinom = 0b1011; // x^3 + x + 1
    int form_polinom;
    printf("Введите образующий полином в 10-тичной форме.\n");
    printf("Пример: x^3+x+1 => 1011\n");
    printf(">>");
    //scanf("%d", &form_polinom);
    form_polinom=10011;
    //form_polinom=1011;
    form_polinom = decimal_to_binary(form_polinom); // Преобразование 1011 => 11
    printf("\n"); // Отступ

    // Количество информационных символов
    int k;
    printf("Введите количество информационных символов.\n");
    printf(">>");
    //scanf("%d", &k);
    //k=3;
    k=9;
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


    //int data[]={12,2,8,1,4,5,10,10,4}; // мой вариант
    int data[]={10,13,5,8,9,2,7,7,7}; // пример от Гали
   // int data[]={13,9,7,14,5,3,2,8,14}; // хуйня со студфайлов
   // int data[]={4,-1,3}; // полином из методички Владимирова
    polinom_t *info_poly=polinom_init(k,gf);
    for(int i=0;i<k;i++)
    {
       polinom_set(info_poly, i, gf_get(gf,data[i]+1));
    }

    polinom_t *encoded=polinom_init(n,gf);
    for(int i=r;i<n;i++)
    {
       polinom_set(encoded, i, info_poly->data[i-r]);
    }

    int temp_array[n]; // массив для хранения остатков от деления
    for (int i=0;i<n;i++) temp_array[i]=encoded->data[i];

    int index=n-1; // контролируем старшую степень делимого
    int value=temp_array[index];
    int round=0; // реализует сдвиг в  temp_array
    int temp_devision[gen_polinom->degree];
    while(index>=gen_polinom->degree-1) // до тех пор пока степень делимого не будет равна степени генераторного полинова, включая её
    {  
        for(int i=0;i<gen_polinom->degree-1;i++)
        {
            // умножаем ген. полином на элемент при старшей степени и складываем его с делимым/остатком
            temp_devision[i]=gf_add(gf_mult(gen_polinom->data[i],value,gf),temp_array[i+k-1-round]);
            temp_array[n-gen_polinom->degree-round+i]=temp_devision[i];
        }//
        round++;
        index--;
        value=temp_devision[gen_polinom->degree-2];
      // index=0;
    }

    for(int i=0;i<r;i++) encoded->data[i]=temp_array[i];
    for(int i=0;i < n;i++) printf("encoded%d = %d\n",i,gf->rev_table[encoded->data[i]]-1);
    printf("\n"); // Отступ

    for (int i=0;i<n;i++) temp_array[i]=encoded->data[i];
    
    //temp_array[0]=gf_get(gf,2+1); // Раскоментировать, чтобы получился полином Владимирова. Результат проверен, сходится с листочком.
    //temp_array[1]=gf_get(gf,2+1);
    //temp_array[2]=gf_get(gf,1+1);
    //temp_array[3]=gf_get(gf,4+1);
    
    index=n-1; // контролируем старшую степень делимого
    
    value=temp_array[index];
    index=n-1;
    round=0;
    while(index>=gen_polinom->degree-1) // до тех пор пока степень делимого не будет равна степени генераторного полинова, включая её
    {  
        for(int i=0;i<gen_polinom->degree-1;i++)
        {
            // умножаем ген. полином на элемент при старшей степени и складываем его с делимым/остатком
            temp_devision[i]=gf_add(gf_mult(gen_polinom->data[i],value,gf),temp_array[i+k-1-round]);
            temp_array[n-gen_polinom->degree-round+i]=temp_devision[i];
        }//
        round++;
        index--;
        value=temp_devision[gen_polinom->degree-2];
      // index=0;
    }


    for(int i=0;i<r;i++) printf("CHECK%d = %d\n",i,gf->rev_table[temp_array[i]]-1);

    // Освобождение памяти, выделенной под данные
    polinom_free(gen_polinom);
    gf_free(gf);

    return 0;
}