
#include <stdio.h>
#include "gflib/gflib.h"
#include "gflib/polinom.h"
#include "gflib/gen_polinom.h"

int count_bits(int num) {
    int count = 0;
    while (num) {
        num >>= 1;
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
    scanf("%d", &form_polinom);
    form_polinom = decimal_to_binary(form_polinom); // Преобразование 1011 => 11
    // printf("\n"); // Отступ

    // Расчет параметров поля
    int m = count_bits(form_polinom)-1; // GF(q^m)

    int n = (1 << m)-1; // Количество кодовых символов
    // Подрязумевается, что n+1 == 2^m

    // Количество информационных символов
    int k;
    while (1) {
        printf("Введите количество информационных символов.\nОно должно находиться в пределах от 1 до %d.\n", n-2);
        printf(">>");
        scanf("%d", &k);
        if (1 <= k && k <= n-2) break;
        printf("Неверный ввод.\n");
    }

    // Остальные параметры

    int r = n - k; // Количество исправляющих символов
    int t = r / 2; // Количество возможных для исправления символов 

    // Формирование поля Галуа
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("gf build error: переданный полином не является неприводимым\n");
        return 1;
    }

    gf_elem_t coded[n];
    int read_data = 1;
    // Чтение данных и проверка на соответствие полю
    while (read_data) {
        read_data = 0;
        printf("Введите полученные кодовые символы (в виде элементов поля) через пробел.\nТаких символов должно быть %d.\n", n);
        printf(">>");
        for (int i = 0; i < n; i++) {
            scanf("%d", &coded[i]);
            if (coded[i] >= gf->total_quantity) read_data = 1;
        }

        if (!read_data) break;
        printf("Получен неверный элемент, повторите попытку.\n");
    }

    polinom_t *encoded = polinom_init(gf, 1);
    polinom_append(encoded, coded, n);

    // Считаем синдромы от b до b+2t-1
    gf_elem_t syndroms[b+2*t-1];
    syndroms[0] = 0;
    for (int i = b; i <= b+2*t-1; i++) {
        syndroms[i] = polinom_call(encoded, gf_get_by_degree(gf, i));
        printf("syndrom %d: %d\n", i, syndroms[i]);
    }


    gf_free(gf);
    return 0;
}