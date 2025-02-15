

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


    int data[]={4,-1,3};

    polinom_t *info_poly=polinom_init(k,gf);
    for(int i=0;i<k;i++)
    {
       polinom_set(info_poly, i, gf_get(gf,data[i]+1));
    }

    polinom_t *encoded=polinom_init(n,gf);
    for(int i=k+1;i<n;i++)
    {
       polinom_set(encoded, i, info_poly->data[i-k-1]);
    }

    printf("degree=%d\n\n",gen_polinom->degree);

    int temp_array[n]; // массив для хранения остатков от деления
    for (int i=0;i<n;i++) temp_array[i]=encoded->data[i];

    int index=n-1; // контролируем старшую степень делимого
    int value=temp_array[index];

    while(index>=gen_polinom->degree-1) // до тех пор пока степень делимого не будет равна степени генераторного полинова, включая её
    {  
        int temp_devision[index];// результат деления. Каждую итерацию размер уменьшается на 1.
        //for (int i=0;i<index;i++) temp_devision[i]=-1;
        // for(int i=index-1;i>=0;i--)
        for(int i=0;i<index;i++)
        {
            // умножаем ген. полином на элемент при старшей степени и складываем его с делимым/остатком
            temp_devision[i]=gf_add(gf_mult(gen_polinom->data[i-(index-gen_polinom->degree+1)],value,gf),temp_array[i]);
            // формируем новый дилитель
            temp_array[i]=temp_devision[i];
        }
        //for(int i=0;i<index;i++) temp_array[i]=temp_devision[i];
        index--;
        value=temp_array[index];
    }

    for(int i=0;i<r;i++) encoded->data[i]=temp_array[i];
    for(int i=0;i < n;i++) printf("encoded%d = %d\n",i,gf->rev_table[encoded->data[i]]-1);
    printf("\n"); // Отступ



    // Освобождение памяти, выделенной под данные
    polinom_free(gen_polinom);
    gf_free(gf);

    return 0;
}