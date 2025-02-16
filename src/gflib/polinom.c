
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gflib.h"

typedef struct polinom polinom_t;

struct polinom
{
    gf_t *gf;           // поле Галуа (на котором основан)
    int capacity;       // емкость полинома (количество ячеек памяти для хранения коэффициентов)
    int degree;         // степень полинома
    gf_elem_t *data;    // массив коэффициентов полинома
};


/*
    Инициализация полинома.
    Можно указать начальную емкость.
*/

polinom_t* polinom_init(gf_t* gf, int capacity) {
    polinom_t *polinom = malloc(sizeof(polinom_t));
    polinom->gf = gf;
    polinom->capacity = capacity;
    polinom->degree = 0;
    polinom->data = calloc(capacity, sizeof(gf_elem_t));
    for (int i = 0; i < capacity; i++) polinom->data[i] = 0;
    return polinom;
}


/*
    Метод для освобождения памяти,
    определенной под полином
*/

void polinom_free(polinom_t *polinom) {
    free(polinom->data);
    free(polinom);
}


/*
    Расширение полинома.
    Предполагается, что полином можно
    будет расширить в связи с ограниченностью
    выделенной памяти.
    По умолчанию размер полинома всегда
    увеличивается в двое.
*/

void polinom_extencion(polinom_t *polinom, int size) {
    gf_elem_t* new_data = calloc(polinom->capacity + size, sizeof(gf_elem_t));
    memcpy(new_data, polinom->data, polinom->capacity * sizeof(gf_elem_t));
    free(polinom->data);
    polinom->data = new_data;
    polinom->capacity += size;
}


/*
    Метод автоматического подсчета степени полинома.
*/

void polinom_calc_degree(polinom_t *polinom) {
    polinom->degree = 0;
    for (int i = 0; i < polinom->capacity; i++) {
        if (polinom->data[i] != 0) polinom->degree = i+1;
    }
}


/*
    Установка значения коэффициента полинома
    по определенному порядковому номеру.
    Реализации для g_elem_t и g_inner_t.
*/

void polinom_set(polinom_t *polinom, int index, gf_elem_t elem) {

    if (index >= polinom->capacity) {
        printf("set error: index is out of capaciry\n");
        printf("index: %d, capacity: %d\n", index, polinom->capacity);
        index = index % polinom->capacity;
    }

    polinom->data[index] = elem;
    if (polinom->degree <= index) polinom->degree = index+1;
    if (elem == 0 && index < polinom->degree) polinom_calc_degree(polinom);
}


/*
    Метод для перемножения двух полиномов
    в установленном поле Галуа.
*/

void polinom_mult(polinom_t *polinom1, polinom_t *polinom2) {
    polinom_t *polinom = polinom_init(polinom1->gf, polinom1->degree + polinom2->degree);
    // TODO: заменить capacity на degree

    if (polinom1->gf != polinom2->gf) {
        printf("polinom mult error: polinoms are not in same GF\n");
        return;
    }

    for (int i = 0; i < polinom1->degree; i++) {
        for (int j = 0; j < polinom2->degree; j++) {
            gf_elem_t mult = gf_mult(polinom1->gf, polinom1->data[i], polinom2->data[j]);
            polinom_set(polinom, i+j, gf_add(polinom->data[i+j], mult));
        }
    }

    if (polinom->degree > polinom1->degree) polinom_extencion(polinom1, polinom->degree - polinom1->degree);

    for (int i = 0; i < polinom->degree; i++) {
        polinom_set(polinom1, i, polinom->data[i]);
    }

    polinom_free(polinom);
}


/*
    Метод для копирования полинома
*/

polinom_t* polinom_copy(polinom_t *polinom1) {
    polinom_t *polinom = polinom_init(polinom1->gf, polinom1->capacity);
    for (int i = 0; i < polinom1->degree; i++) {
        polinom_set(polinom, i, polinom1->data[i]);
    }
    return polinom;
}


/*
    Метод для вывода коэффициентов полинома
    через функцию printf
*/

void polinom_print(polinom_t *polinom) {
    int is_first = 1;
    for (int i = 0; i < polinom->degree; i++) {
        if (polinom->data[i] == 0) continue;
        if (!is_first) printf(" + ");
        else is_first = 0;
        printf("e^%d*x^%d", polinom->gf->rev_table[polinom->data[i]]-1, i);
    }
    if (is_first) printf("0");
    printf("\n");
}