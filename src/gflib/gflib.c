
#include <stdio.h>
#include <stdlib.h>

/*
    GF(q^m), где q - простое число, m - степень поля
    В данной реализации используется направление
    полиномов от большей к меньше степени,
    то есть a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x^1 + f 
*/

typedef struct gf gf_t;
typedef int gf_inner_t;
typedef struct gf_elem gf_elem_t;

struct gf
{
    int characteristic;     // q
    int power;              // m
    int total_quantity;     // q^m
    gf_inner_t *table;       // таблица элементов поля
    int forming_polinom;    // полином для построения
    int mask;               // маска для mod операции
};

struct gf_elem
{
    gf_inner_t value;
    gf_t* gf;
};


/*
    Инициализация поля Галуа
    В текущей реализации q = 2
    Для большей универсальности необходима
    реализация возвредения в целочисленную
    степень.
*/

gf_t* gf_init(int power) {
    gf_t *gf = malloc(sizeof(gf_t));
    gf->characteristic = 2;
    gf->power = power;
    gf->total_quantity = 1 << power;
    gf->table = 0;
    return gf;
}


/*
    Построение таблицы элементов поля
    Пример для q = 2, m = 3, forming_polinom = 1011 (x^3 + x + 1)

    Получается следующая таблица:
    -----------------------------
    | id | primitive | bits | nums |
    | 0  |    0      | 000  |  0   |
    | 1  |    1      | 001  |  1   |
    | 2  |    e      | 010  |  2   |
    | 3  |   e^2     | 100  |  4   |
    | 4  |   e+1     | 011  |  3   |
    | 5  |  e^2+e    | 110  |  6   |
    | 6  | e^2+e+1   | 111  |  7   |
    | 7  |  e^2+1    | 101  |  5   |
*/

void gf_build(gf_t *gf, int polinom) {
    gf->table = calloc(gf->total_quantity, sizeof(int));
    gf->forming_polinom = polinom;

    gf->mask = 1 << gf->power;
    gf->table[0] = 0;
    gf->table[1] = 1;

    for (int i = 2; i < gf->total_quantity; i++) {
        gf->table[i] = gf->table[i-1] << 1;
        if(gf->table[i] & gf->mask) {
            gf->table[i] ^= polinom;
        }
    }
}


/*
    Получение элемента поля
*/

gf_elem_t gf_get(gf_t *gf, int id) {
    gf_elem_t tmp = { gf->table[id], gf };
    return tmp;
}


/*
    Сложение элементов поля
*/

gf_elem_t gf_add(gf_elem_t a, gf_elem_t b) {

    gf_elem_t tmp = { a.value ^ b.value, a.gf };

    if (a.gf != b.gf) {
        tmp.value = 0;
        tmp.gf = 0;
        printf("add error: GF a is not in GF b\n");
    }

    return tmp;
}


/*
    Умножение элементов поля
*/

gf_elem_t gf_mult(gf_elem_t a, gf_elem_t b) {
    gf_elem_t res = { 0, a.gf };

    if (a.gf != b.gf) {
        res.gf = 0;
        printf("mult error: GF a is not in GF b\n");
    }

    for (int i = 0; i < a.gf->power; i++) {
        if (b.value & (1 << i)) {
            res.value ^= a.value << i;
            if (res.value & a.gf->mask) {
                res.value ^= a.gf->forming_polinom;
            }
        }
    }
    return res;
}

/*
    Печать элементов поля
*/

void gf_print(gf_t *gf) {
    for (int i = 0; i < gf->total_quantity; i++) {
        printf("%d: %d\n", i, gf->table[i]);
    }
}


/*
    Освобождение поля
*/

void gf_free(gf_t *gf) {
    free(gf);
}