
#include <stdio.h>
#include <stdlib.h>

/*
    GF(q^m), где q - простое число, m - степень поля
    В данной реализации используется направление
    полиномов от большей к меньше степени,
    то есть a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x^1 + f 
*/

struct gf
{
    int characteristic;     // q
    int power;              // m
    int total_quantity;     // q^m
    int *table;             // таблица элементов поля
    int polinom;            // полином для построения
};

typedef struct gf gf_t;
typedef int gf_elem_t;

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
    Пример для q = 2, m = 3, polinom = 1011 (x^3 + x + 1)

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
    gf->polinom = polinom;

    int mask = 1 << gf->power;
    gf->table[0] = 0;
    gf->table[1] = 1;

    for (int i = 2; i < gf->total_quantity; i++) {
        gf->table[i] = gf->table[i-1] << 1;
        if(gf->table[i] & mask) {
            gf->table[i] ^= polinom;
        }
    }
}

void gf_print(gf_t *gf) {
    for (int i = 0; i < gf->total_quantity; i++) {
        printf("%d: %d\n", i, gf->table[i]);
    }
}

void gf_free(gf_t *gf) {
    free(gf);
}