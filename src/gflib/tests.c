
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
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    int data[] = { 0b000, 0b001, 0b010, 0b100, 0b011, 0b110, 0b111, 0b101 };
    int length = sizeof(data) / sizeof(data[0]);
    for (int i = 0; i < length; i++) {
        if (data[i] != gf->table[i]) {
            printf("\t\tgf build error: расхождение при построении поля Галуа\n" \
                "\t\treal[%d] == %d\n" \
                "\t\tgf->table[%d] == %d\n",
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
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    gf_elem_t a = gf_get_by_degree(gf, 4);
    gf_elem_t b = gf_get_by_id(gf, 5);
    gf_elem_t c = gf_get_by_id(gf, 0);

    if (a != b) {
        printf("\t\tgf get error: %d != %d\n", a, b);
        return 1;
    }

    if (c != 0) {
        printf("\t\tgf get error: %d != 0\n", c);
        return 1;
    }

    gf_free(gf);
    return 0;
}

int test4() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    // Сложение, умножение и возведение в степень
    gf_elem_t a = gf_get_by_degree(gf, 4);
    gf_elem_t b = gf_get_by_degree(gf, 5);
    gf_elem_t c = gf_add(a, b);
    gf_elem_t d = gf_mult(gf, a, b);
    gf_elem_t e = gf_pow(gf, a, 20);

    if (c != gf_get_by_degree(gf, 0)) {
        printf("\t\tgf add error: %d + %d != %d\n", a, b, c);
        return 1;
    }

    if (d != gf_get_by_degree(gf, 9)) {
        printf("\t\tgf mult error: %d * %d != %d\n", a, b, d);
        return 1;
    }

     if (e != gf_get_by_degree(gf, 3)) {
        printf("\t\tgf pow error: %d^%d != %d\n", a, 20, e);
        return 1;
    }

    // Деление
    a = gf_get_by_degree(gf, 0);
    b = gf_get_by_degree(gf, 3);
    c = gf_get_by_degree(gf, 5);
    d = gf_div(gf, a, b);
    e = gf_div(gf, c, b);

    if (d != gf_get_by_degree(gf, 3)) {
        printf("\t\tgf div error (no carry): %d / %d != %d\n", a, b, d);
        return 1;
    }

    if (e != gf_get_by_degree(gf, 5)) {
        printf("\t\tgf div error (with carry): %d / %d != %d\n", b, c, e);
        return 1;
    }


    gf_free(gf);
    return 0;
}

int test5() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    // Инициализация полинома
    int size = 10;
    polinom_t *polinom = polinom_init(gf, size);
    polinom->data[size-1] = 1;
    polinom_extencion(polinom, size*3);

    if (polinom->capacity != size*4) {
        printf("\t\tpolinom extension error: capacity != %d\n", size*4);
        return 1;
    }

    int i = 0;
    while (i < polinom->capacity && polinom->data[i] == 0) i++;

    if (i != size-1) {
        printf("\t\tpolinom extension error: did not found '1' in data\n");
        return 1;
    }

    i++;
    while (i < polinom->capacity && polinom->data[i] == 0) i++;

    if (i != size*4) {
        printf("\t\tpolinom extension error: did not found end of capacity\n");
        return 1;
    }

    polinom_free(polinom);
    gf_free(gf);
    return 0;
}

int test6() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    // Инициализация полинома
    polinom_t *polinom = polinom_init(gf, 10);

    if (polinom->degree != 0) {
        printf("\t\tpolinom degree error: degree != 0\n");
        return 1;
    }

    polinom->data[1] = 1;
    polinom_calc_degree(polinom);

    if (polinom->degree != 2) {
        printf("\t\tpolinom degree error: degree != 2\n");
        return 1;
    }

    polinom->data[1] = 0;
    polinom->data[2] = 1;
    polinom->data[5] = 1;
    polinom_calc_degree(polinom);

    if (polinom->degree != 6) {
        printf("\t\tpolinom degree error: degree != 6\n");
        return 1;
    }

    polinom->data[5] = 0;
    polinom_calc_degree(polinom);

    if (polinom->degree != 3) {
        printf("\t\tpolinom degree error: degree != 3\n");
        return 1;
    }

    polinom_free(polinom);
    gf_free(gf);
    return 0;
}

int test7() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    // Инициализация полинома
    polinom_t *polinom = polinom_init(gf, 10);

    if (polinom->degree != 0) {
        printf("\t\tpolinom degree error: degree != 0\n");
        return 1;
    }

    polinom_set(polinom, 1, 1);

    if (polinom->degree != 2) {
        printf("\t\tpolinom set error: degree != 2\n");
        return 1;
    }

    polinom_set(polinom, 1, 0);
    polinom_set(polinom, 2, 1);
    polinom_set(polinom, 5, 1);

    if (polinom->degree != 6) {
        printf("\t\tpolinom set error: degree != 6\n");
        return 1;
    }

    polinom_set(polinom, 5, 0);
    if (polinom->degree != 3) {
        printf("\t\tpolinom set error: degree != 3\n");
        return 1;
    }

    polinom_free(polinom);
    gf_free(gf);
    return 0;
}

int test8() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    // Инициализация полинома
    polinom_t *polinom = polinom_init(gf, 1);
    gf_elem_t data[] = {
        gf_get_by_degree(gf, 0),
        gf_get_by_degree(gf, 1),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 4),
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 6),
    };
    polinom_append(polinom, data, sizeof(data)/sizeof(data[0]));
    
    if (polinom->degree != 7 && polinom->capacity != 7) {
        printf("\t\tpolinom append error: polinom did not extended correctly\n");
        return 1;
    }

    for (int i = 0; i < polinom->degree; i++) {
        if (data[i] != polinom->data[i]) {
            printf("\t\tpolinom append error: append does not correct (%d != %d)\n", data[i], polinom->data[i]);
            return 1;
        }
    }

    polinom_free(polinom);
    gf_free(gf);
    return 0;
}

int test9() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    polinom_t *polinom1 = polinom_init(gf, 1);
    gf_elem_t data1[] = {
        gf_get_by_degree(gf, 0),
        gf_get_by_degree(gf, 1),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 4),
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 6),
    };
    polinom_append(polinom1, data1, sizeof(data1)/sizeof(data1[0]));

    polinom_t *polinom2 = polinom_init(gf, 1);
    gf_elem_t data2[] = {
        gf_get_by_degree(gf, 6),
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 4),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 1),
        gf_get_by_degree(gf, 0),
        gf_get_by_degree(gf, 3)
    };
    polinom_append(polinom2, data2, sizeof(data2)/sizeof(data2[0]));

    polinom_add(polinom1, polinom2);
    polinom_free(polinom2);
    
    gf_elem_t data3[] = {
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 6),
        gf_get_by_degree(gf, 1),
        0,
        gf_get_by_degree(gf, 1),
        gf_get_by_degree(gf, 6),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
    };

    if (polinom1->degree != 8 && polinom1->capacity != 8) {
        printf("\t\tpolinom add error: polinom did not extended correctly\n");
        return 1;
    }

    for (int i = 0; i < polinom1->degree; i++) {
        if (data3[i] != polinom1->data[i]) {
            printf("\t\tpolinom add error: append does not correct (%d != %d)\n", data3[i], polinom1->data[i]);
            return 1;
        }
    }
    

    polinom_free(polinom1);
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

    // Тест 3: получение элементов поля Галуа
    printf("\ttest 3: получение элементов поля Галуа\n");
    if (test3()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 4: операции с элементами поля Галуа
    printf("\ttest 4: операции с элементами поля Галуа\n");
    if (test4()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    printf("\n");

    printf("Тестирование polinom\n");

    // Тест 5: проверка полинома на корректную инициализацию и расширение
    printf("\ttest 5: проверка полинома на корректную инициализацию и расширение\n");
    if (test5()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 6: проверка функции отсчета степени полинома
    printf("\ttest 6: проверка функции отсчета степени полинома\n");
    if (test6()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 7: проверка фукнции polinom_set
    printf("\ttest 7: проверка фукнции polinom_set\n");
    if (test7()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 8: проверка фукнции polinom_append
    printf("\ttest 8: проверка фукнции polinom_append\n");
    if (test8()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 9: проверка фукнции polinom_add
    printf("\ttest 9: проверка фукнции polinom_add\n");
    if (test9()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    return 0;
}