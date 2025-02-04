

#include <stdio.h>
#include "gflib/gflib.h"
#include "gflib/polinom.h"
#include "gflib/gen_polinom.h"

int main() {

    // Подготовка к использованию
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1

    gf_t *gf = gf_init(m);
    gf_build(gf, form_polinom);
    gf_print(gf); // Отобразим элементы поля Галуа


    // Пример 1: получаем элемента поля Галуа
    // и производим сложение и умножение
    gf_elem_t a = gf_get(gf, 4);
    gf_elem_t b = gf_get(gf, 5);
    gf_elem_t c = gf_add(a, b);
    gf_elem_t d = gf_mult(a, b);

    printf("%d + %d = %d\n", a.value, b.value, c.value);
    printf("%d * %d = %d\n", a.value, b.value, d.value);


    // Пример 2: Умножение полиномов
    // Полином 1: e^1 + x
    polinom_t *polinom1 = polinom_init(5, gf);
    polinom_set(polinom1, 0, gf_get(gf, 2));
    polinom_set(polinom1, 1, gf_get(gf, 1));
    polinom_print(polinom1);
    printf("degree: %d\n", polinom1->degree);
    printf("capacity: %d\n", polinom1->capacity);

    // Полином 2: e^2 + x
    polinom_t *polinom2 = polinom_init(5, gf);
    polinom_set(polinom2, 0, gf_get(gf, 3));
    polinom_set(polinom2, 1, gf_get(gf, 1));
    polinom_print(polinom2);
    printf("degree: %d\n", polinom1->degree);
    printf("capacity: %d\n", polinom1->capacity);

    // Полином 3: результат умножения П1 на П2
    polinom_t *polinom3 = polinom_mult(polinom1, polinom2);
    polinom_print(polinom3);
    printf("degree: %d\n", polinom3->degree);
    printf("capacity: %d\n", polinom3->capacity);


    // Пример 3: формирование порождающего полинома
    polinom_t *gen_polinom = generating_polinom(gf, 1, 2);
    polinom_print(gen_polinom);


    // Освобождение памяти, выделенной под данные
    polinom_free(polinom1);
    polinom_free(polinom2);
    polinom_free(polinom3);
    polinom_free(gen_polinom);
    gf_free(gf);

    return 0;
}