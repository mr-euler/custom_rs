

#include <stdio.h>
#include "gflib/gflib.h"
#include "gflib/polinom.h"
#include "gflib/gen_polinom.h"

int main() {

    // Подготовка к использованию
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1

    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("gf build error: polinom is not primitive\n");
        return 1;
    }
    gf_print_bin(gf); // Отобразим элементы поля Галуа
    printf("\n");


    // Пример 1: получаем элемента поля Галуа
    // и производим сложение и умножение
    gf_elem_t a = gf_get_by_id(gf, 4);
    gf_elem_t b = gf_get_by_id(gf, 5);
    gf_elem_t c = gf_add(a, b);
    gf_elem_t d = gf_mult(gf, a, b);

    printf("%d + %d = %d\n", a, b, c);
    printf("%d * %d = %d\n", a, b, d);
    printf("\n");


    // Пример 2: Умножение полиномов
    // Полином 1: e^1 + x
    polinom_t *polinom1 = polinom_init(gf, 5);
    polinom_set(polinom1, 0, gf_get_by_id(gf, 2));
    polinom_set(polinom1, 1, gf_get_by_id(gf, 1));
    polinom_print(polinom1);
    printf("degree: %d\n", polinom1->degree);
    printf("capacity: %d\n", polinom1->capacity);
    printf("\n");

    // Полином 2: e^2 + x
    polinom_t *polinom2 = polinom_init(gf, 5);
    polinom_set(polinom2, 0, gf_get_by_id(gf, 3));
    polinom_set(polinom2, 1, gf_get_by_id(gf, 1));
    polinom_print(polinom2);
    printf("degree: %d\n", polinom1->degree);
    printf("capacity: %d\n", polinom1->capacity);
    printf("\n");

    // Полином 3: результат умножения П1 на П2
    polinom_mult(polinom1, polinom2);
    polinom_print(polinom1);
    printf("degree: %d\n", polinom1->degree);
    printf("capacity: %d\n", polinom1->capacity);
    printf("\n");


    // Пример 3: формирование порождающего полинома
    polinom_t *gen_polinom = generating_polinom(gf, 1, 2);
    printf("Генераторный полином:\n");
    polinom_print(gen_polinom);
    printf("\n");

    // Освобождение памяти, выделенной под данные
    polinom_free(polinom1);
    polinom_free(polinom2);
    polinom_free(gen_polinom);
    gf_free(gf);

    return 0;
}