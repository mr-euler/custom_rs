
#include <stdio.h>
#include "gflib.h"
#include "polinom.h"
#include "gen_polinom.h"

int test1() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1001; // x^3 + 1
    gf_t *gf = gf_init(m);
    if(!gf_build(gf, form_polinom)) {
        gf_free(gf);
        return 1;
    }
    return 0;
}

int test2() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("gf build error: полином не является неприводимый\n");
        return 1;
    }

    int data[] = { 0b000, 0b001, 0b010, 0b100, 0b011, 0b110, 0b111, 0b101 };
    int length = sizeof(data) / sizeof(data[0]);
    for (int i = 0; i < length; i++) {
        if (data[i] != gf->table[i]) {
            printf("gf build error: расхождение при построении поля Галуа\n" \
                "real[%d] == %d\n" \
                "gf->table[%d] == %d\n",
                i, data[i], i, gf->table[i]);
            return 1;
        }
    }
    gf_free(gf);
    return 0;
}

int test3() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("gf build error: полином не является неприводимый\n");
        return 1;
    }

    gf_elem_t a = gf_get_by_degree(gf, 4);
    gf_elem_t b = gf_get_by_degree(gf, 5);
    gf_elem_t c = gf_add(a, b);
    gf_elem_t d = gf_mult(gf, a, b);

    if (c != gf_get_by_degree(gf, 0)) {
        printf("gf add error: %d + %d != %d\n", a, b, c);
        return 1;
    }

    if (d != gf_get_by_degree(gf, 9)) {
        printf("gf mult error: %d * %d != %d\n", a, b, c);
        return 1;
    }

    gf_free(gf);
    return 0;
}

int main() {

    printf("Тестирование gflib\n");

    // Тест 1: проверка полинома на приводимость
    printf("\ttest 1: проверка полинома на приводимость\n");
    if (test1()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 2: построение поля Галуа
    printf("\ttest 2: построение поля Галуа\n");
    if (test2()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 3: операции с элементами поля Галуа
    printf("\ttest 3: операции с элементами поля Галуа\n");
    if (test3()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    printf("\n");

    return 0;
}