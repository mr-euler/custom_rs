

#include <stdio.h>
#include "gflib/gflib.h"
#include "gflib/polinom.h"
#include "gflib/gen_polinom.h"

int main() {

    // Параметры кода РС

    // Константы

    int n = 7; // Количество кодовых символов
    // Подрязумевается, что n+1 == 2^m

    int b = 1; // Параметр кодирования

    int m = 0; // GF(q^m)
    { // Вычисление степени поля по количеству кодовых символов
        int tmp = n;
        while (tmp) {
            tmp >>= 1;
            m++;
        }

    }

    // Динамические параметры

    // Полином
    // int form_polinom = 0b1011; // x^3 + x + 1
    int form_polinom;
    printf("Введите образующий полином в 10-тичной форме.\n");
    printf("Пример: x^3+x+1 => 0b1011 => 11\n");
    printf(">>");
    scanf("%d", &form_polinom);
    printf("\n"); // Отступ

    // Количество информационных символов
    int k;
    printf("Введите количество информационных символов.\n");
    printf(">>");
    scanf("%d", &k);
    printf("\n"); // Отступ

    // Остальные параметры

    int r = n - k; // Количество исправляющих символов
    int t = r / 2; // Количество возможных для исправления символов 

    // Формирование поля Галуа
    gf_t *gf = gf_init(m);
    gf_build(gf, form_polinom);
    // gf_print(gf); // Отобразим элементы поля Галуа


    // Формирование порождающего полинома
    polinom_t *gen_polinom = generating_polinom(gf, b, t);
    printf("Порождающий полином:\n");
    polinom_print(gen_polinom);


    // Освобождение памяти, выделенной под данные
    polinom_free(gen_polinom);
    gf_free(gf);

    return 0;
}