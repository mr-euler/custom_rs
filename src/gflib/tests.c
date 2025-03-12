
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
            printf("\t\tpolinom add error: add does not correct (%d != %d)\n", data3[i], polinom1->data[i]);
            return 1;
        }
    }
    

    polinom_free(polinom1);
    gf_free(gf);
    return 0;
}

int test10() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    polinom_t *polinom1 = polinom_init(gf, 1);
    gf_elem_t data1[] = {
        gf_get_by_degree(gf, 1),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
    };
    polinom_append(polinom1, data1, sizeof(data1)/sizeof(data1[0]));

    polinom_t *polinom2 = polinom_init(gf, 1);
    gf_elem_t data2[] = {
        gf_get_by_degree(gf, 4),
        gf_get_by_degree(gf, 6),
    };
    polinom_append(polinom2, data2, sizeof(data2)/sizeof(data2[0]));

    polinom_mult(polinom1, polinom2);
    polinom_free(polinom2);
    
    gf_elem_t data3[] = {
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 2),
    };

    if (polinom1->degree != 4 && polinom1->capacity != 4) {
        printf("\t\tpolinom mult error: polinom did not extended correctly\n");
        return 1;
    }

    for (int i = 0; i < polinom1->degree; i++) {
        if (data3[i] != polinom1->data[i]) {
            printf("\t\tpolinom mult error: mult does not correct (%d != %d)\n", data3[i], polinom1->data[i]);
            return 1;
        }
    }
    

    polinom_free(polinom1);
    gf_free(gf);
    return 0;
}

int test11() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    polinom_t *polinom1 = polinom_init(gf, 1);
    gf_elem_t data1[] = {
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 2),
    };
    polinom_append(polinom1, data1, sizeof(data1)/sizeof(data1[0]));

    gf_elem_t result = polinom_call(polinom1, gf_get_by_degree(gf, 2));

    if (result != gf_get_by_degree(gf, 1)) {
        printf("\t\tpolinom call error: wrong answer (%d != %d)\n", result, gf_get_by_degree(gf, 1));
        return 1;
    }
    

    polinom_free(polinom1);
    gf_free(gf);
    return 0;
}

int test12() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    polinom_t *polinom1 = polinom_init(gf, 1);
    gf_elem_t data1[] = {
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 2),
    };
    polinom_append(polinom1, data1, sizeof(data1)/sizeof(data1[0]));

    polinom_t *polinom2 = polinom_copy(polinom1);

    if (polinom1->capacity != polinom2->capacity) {
        printf("\t\tpolinom copy error: wrong capacity\n");
        return 1;
    }

    if (polinom1->degree != polinom2->degree) {
        printf("\t\tpolinom copy error: wrong degree\n");
        return 1;
    }

    for (int i = 0; i < polinom1->degree; i++) {
        if (data1[i] != polinom1->data[i]) {
            printf("\t\tpolinom copy error: wrong element (%d != %d)\n", data1[i], polinom1->data[i]);
        return 1;
        }
    }
    

    polinom_free(polinom1);
    polinom_free(polinom2);
    gf_free(gf);
    return 0;
}

int test13() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    polinom_t *polinom1 = polinom_init(gf, 1);
    gf_elem_t data1[] = {
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 2),
    };
    polinom_append(polinom1, data1, sizeof(data1)/sizeof(data1[0]));

    polinom_clear(polinom1);

    if (polinom1->degree != 0) {
        printf("\t\tpolinom clear error: degree != 0 (%d)\n", polinom1->degree);
        return 1;
    }

    for (int i = 0; i < polinom1->capacity; i++) {
        if (polinom1->data[i] != 0) {
            printf("\t\tpolinom clear error: wrong element (%d != 0)\n", polinom1->data[i]);
            return 1;
        }
    }
    
    polinom_free(polinom1);
    gf_free(gf);
    return 0;
}

int test14() {
    int m = 4; // GF(q^m)
    int form_polinom = 0b10011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    polinom_t *info_polinom = polinom_init(gf, 1);
    // gf_elem_t data1[] = { // Галин
    //     gf_get_by_degree(gf, 10),
    //     gf_get_by_degree(gf, 13),
    //     gf_get_by_degree(gf, 5),
    //     gf_get_by_degree(gf, 8),
    //     gf_get_by_degree(gf, 9),
    //     gf_get_by_degree(gf, 2),
    //     gf_get_by_degree(gf, 7),
    //     gf_get_by_degree(gf, 7),
    //     gf_get_by_degree(gf, 7),
    // };

    // int data[]={13,9,7,14,5,3,2,8,14};
    gf_elem_t data1[] = { // "Хуйня со студфайла" (c) Жак Фрэско
        gf_get_by_degree(gf, 13),
        gf_get_by_degree(gf, 9),
        gf_get_by_degree(gf, 7),
        gf_get_by_degree(gf, 14),
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 8),
        gf_get_by_degree(gf, 14),
    };
    polinom_append(info_polinom, data1, sizeof(data1)/sizeof(data1[0]));
    
    // printf("Информационный полином: ");
    // polinom_print(info_polinom);
    // printf("\n");

    int n = (1 << m)-1; // Количество кодовых символов
    // Подрязумевается, что n+1 == 2^m

    int b = 1;
    int k = 9;
    int r = n - k; // Количество исправляющих символов
    int t = r / 2; // Количество возможных для исправления символов 

    // polinom_t *tmp_polinom = polinom_init(gf, r+1);
    // polinom_set(tmp_polinom, tmp_polinom->capacity-1, 1);

    polinom_t *encoded = polinom_copy(info_polinom);
    polinom_right_shift(encoded, r);
    // polinom_mult(encoded, tmp_polinom);

    // printf("Информационный полином (сдвиг): ");
    // polinom_print(encoded);
    // printf("\n");


    polinom_t *gen_polinom = generating_polinom(gf, b, t);

    polinom_t *mod = polinom_copy(encoded);
    polinom_mod(mod, gen_polinom);

    // printf("Остаток: ");
    // polinom_print(mod);
    // printf("\n");

    polinom_add(encoded, mod);

    // printf("Итоговый кодовый полином: ");
    // polinom_print(encoded);
    // printf("\n");

    polinom_mod(encoded, gen_polinom);

    // printf("Проверка: ");
    // polinom_print(encoded);
    // printf("\n");
    
    for (int i = 0; i < encoded->degree; i++) {
        if (encoded->data[i] != 0) {
            printf("\t\tpolinom mod error: wrong element (%d != 0)\n", encoded->data[i]);
            return 1;
        }
    }

    
    polinom_free(info_polinom);
    // polinom_free(tmp_polinom);
    polinom_free(encoded);
    polinom_free(gen_polinom);
    polinom_free(mod);
    gf_free(gf);
    return 0;
}

int test15() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    polinom_t *polinom1 = polinom_init(gf, 6);
    gf_elem_t data1[] = {
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 2),
    };
    polinom_append(polinom1, data1, sizeof(data1)/sizeof(data1[0]));

    polinom_right_shift(polinom1, 4);

    if (polinom1->degree != 8) {
        printf("\t\tpolinom right shift error: degree != 8 (%d)\n", polinom1->degree);
        return 1;
    }

    gf_elem_t check[] = {
        0,
        0,
        0,
        0,
        gf_get_by_degree(gf, 5),
        gf_get_by_degree(gf, 2),
        gf_get_by_degree(gf, 3),
        gf_get_by_degree(gf, 2),
    };

    for (int i = 0; i < polinom1->degree; i++) {
        if (polinom1->data[i] != check[i]) {
            printf("\t\tpolinom right shift error: wrong element (%d != 0)\n", polinom1->data[i]);
            return 1;
        }
    }
    
    polinom_free(polinom1);
    gf_free(gf);
    return 0;
}

int test4_1() {
    int m = 3; // GF(q^m)
    int form_polinom = 0b1011; // x^3 + x + 1
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("\t\tgf build error: полином не является неприводимый\n");
        return 1;
    }

    gf_elem_t neg = gf_neg(gf, gf_get_by_degree(gf, 1));

    if (neg != gf_get_by_degree(gf, 6)) {
        printf("\t\tgf negative error: e^-1 != e^6 (%d)\n", neg);
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

    // Тест 4.1: операция получения обратного элемента поля Галуа
    printf("\ttest 4.1: операция получения обратного элемента поля Галуа\n");
    if (test4_1()) {
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

    // Тест 10: проверка фукнции polinom_add
    printf("\ttest 10: проверка фукнции polinom_mult\n");
    if (test10()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 11: проверка фукнции polinom_call
    printf("\ttest 11: проверка фукнции polinom_call\n");
    if (test11()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 12: проверка фукнции polinom_copy
    printf("\ttest 12: проверка фукнции polinom_copy\n");
    if (test12()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 13: проверка фукнции polinom_clear
    printf("\ttest 13: проверка фукнции polinom_clear\n");
    if (test13()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 14: проверка фукнции polinom_mod
    printf("\ttest 14: проверка фукнции polinom_mod\n");
    if (test14()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    // Тест 15: проверка фукнции polinom_right_shift
    printf("\ttest 15: проверка фукнции polinom_right_shift\n");
    if (test15()) {
        printf("\t\tне пройден\n");
        return 1;
    }
    printf("\t\tпройден\n");

    return 0;
}