

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
    // printf("\n"); // Отступ

    // Остальные параметры

    int r = n - k; // Количество исправляющих символов
    int t = r / 2; // Количество возможных для исправления символов 

    // Формирование поля Галуа
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("gf build error: переданный полином не является неприводимым\n");
        return 1;
    }
    // gf_print(gf); // Отобразим элементы поля Галуа

    gf_elem_t info[k];
    int read_data = 1;
    // Чтение данных и проверка на соответствие полю
    while (read_data) {
        read_data = 0;
        printf("Введите информационные символы (в виде элементов поля) через пробел.\n");
        printf(">>");
        for (int i = 0; i < k; i++) {
            scanf("%d", &info[i]);
            if (info[i] >= gf->total_quantity) read_data = 1;
        }

        if (!read_data) break;
        printf("Получен неверный элемент, повторите попытку.\n");
    }

    // Преобразовываем информационные символы в полином
    polinom_t *info_polinom = polinom_init(gf, 1);
    polinom_append(info_polinom, info, sizeof(info)/sizeof(info[0]));

    // Формирование порождающего полинома
    polinom_t *gen_polinom = generating_polinom(gf, b, t);

    // Формируем кодовый полином
    polinom_t *encoded = polinom_copy(info_polinom);
    polinom_right_shift(encoded, r); // И сдвигаем его на r вправо

    // Получаем остаток от деления
    polinom_t *mod = polinom_copy(encoded);
    polinom_mod(mod, gen_polinom);
    polinom_add(encoded, mod); // И "склеиваем" его с кодовым полиномом

    // Отображаем кодовый полином
    printf("Итоговый кодовый вектор.\n[");
    int c = 0;
    for (int i = 0; i < encoded->capacity; i++) {
        printf(" %d", encoded->data[i]);
        c++;
    }

    // Дополняем нулями
    while (c < n) {
        printf(" 0");
        c++;
    }

    printf(" ]\n");

    // Освобождение памяти, выделенной под данные
    polinom_free(mod);
    polinom_free(encoded);
    polinom_free(info_polinom);
    polinom_free(gen_polinom);
    gf_free(gf);

    return 0;
}