
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

    if (capacity < 1) {
        printf("polinom init warning: capacity < 1\n");
        capacity = 1;
    }

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
    int old_size = polinom->capacity;
    gf_elem_t* new_data = calloc(polinom->capacity + size, sizeof(gf_elem_t));
    memcpy(new_data, polinom->data, polinom->capacity * sizeof(gf_elem_t));
    free(polinom->data);
    polinom->data = new_data;
    polinom->capacity += size;
    for (int i = old_size; i < polinom->capacity; i++) polinom->data[i] = 0;
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
    Метод для добалвения множества коэффициентов
*/

void polinom_append(polinom_t *polinom, int arr[], int size) {
    if (polinom->capacity < size) polinom_extencion(polinom, size - polinom->capacity);
    for (int i = 0; i < size; i++) {
        if (!arr[i]) continue;
        polinom_set(polinom, i, arr[i]);
    }
}


/*
    Метод для сложение полиномов
    Результат сложения записывается
    в первый переданный полином
*/

void polinom_add(polinom_t *polinom1, polinom_t *polinom2) {
    if (polinom1->gf != polinom2->gf) {
        printf("polinom add error: polinoms are not in same GF\n");
        return;
    }

    int i = 0;
    while (i < polinom1->degree && i < polinom2->degree) {
        polinom_set(polinom1, i, gf_add(polinom1->data[i], polinom2->data[i]));
        i++;
    }

    if (i < polinom1->degree) return;
    if (i < polinom2->degree) {
        if (polinom1->capacity < polinom2->degree)
            polinom_extencion(polinom1, polinom2->degree - polinom1->capacity);
        while (i < polinom2->degree) {
            polinom_set(polinom1, i, polinom2->data[i]);
            i++;
        }
    }
}


/*
    Метод для перемножения двух полиномов
    в установленном поле Галуа.
*/

void polinom_mult(polinom_t *polinom1, polinom_t *polinom2) {
    polinom_t *polinom = polinom_init(polinom1->gf, polinom1->degree + polinom2->degree - 1);
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

    polinom1->capacity = polinom->capacity;
    polinom1->degree = polinom->degree;
    free(polinom1->data);
    polinom1->data = polinom->data;
    free(polinom);
}


/*
    Метод для "вызова" полинома
    То есть, берется коэффициент
    (аргумент gf_elem_t (e^m))
    и подставляется за место
    каждого X^n, вычисляется (e^m)^n
    и суммируется, а по завершени
    возвращается

*/

gf_elem_t polinom_call(polinom_t *polinom, gf_elem_t elem) {
    gf_elem_t result = 0;
    for (int i = 0; i < polinom->degree; i++) {
        result = gf_add(
            result,
            gf_mult(
                polinom->gf,
                polinom->data[i],
                gf_pow(
                    polinom->gf,
                    elem,
                    i
                )
            )
        );
    }
    return result;
}


/*
    Метод для копирования полинома
*/

polinom_t* polinom_copy(polinom_t *polinom1) {
    polinom_t *polinom2 = polinom_init(polinom1->gf, polinom1->capacity);
    for (int i = 0; i < polinom1->degree; i++) {
        if (polinom1->data[i])
            polinom_set(polinom2, i, polinom1->data[i]);
    }
    return polinom2;
}


/*
    Метод для очистки полинома
*/

void polinom_clear(polinom_t *polinom) {
    for (int i = 0; i < polinom->capacity; i++) {
        polinom->data[i] = 0;
    }
    polinom->degree = 0;
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


/*
    Метод для получения остатка от
    деления полинома на полином
*/

void polinom_mod(polinom_t *dividend, polinom_t *divisor) {
    while (dividend->degree >= divisor->degree) {
        polinom_t *tmp = polinom_copy(divisor);
        polinom_t *shift = polinom_init(dividend->gf, dividend->degree - divisor->degree + 1);
        gf_elem_t multiplier = gf_div(
            divisor->gf,
            divisor->data[divisor->degree-1],
            dividend->data[dividend->degree-1]
        );

        // gf_elem_t check = gf_mult(dividend->gf, multiplier, divisor->data[divisor->degree-1]);
        // if (check != dividend->data[dividend->degree-1]) {
        //     printf("WA ");
        //     gf_elem_print(divisor->gf, divisor->data[divisor->degree-1]); printf(" ");
        //     gf_elem_print(dividend->gf, dividend->data[dividend->degree-1]); printf(" ");
        //     gf_elem_print(divisor->gf, multiplier); printf("\n");
        //     return;
        // }

        // Данные для сохранения целой части от деления:
        // коэффициент: multiplier
        // степень: dividend->degree - divisor->degree

        polinom_set(shift, shift->capacity-1, multiplier);
        polinom_mult(tmp, shift);
        polinom_add(dividend, tmp);
        polinom_free(tmp);
        polinom_free(shift);
    }
}

