
#include <stdio.h>
#include <stdlib.h>

/*
    GF(q^m), где q - простое число, m - степень поля.
    В данной реализации используется направление
    полиномов от большей к меньше степени,
    то есть a*e^5 + b*e^4 + c*e^3 + d*e^2 + f*e^1 + g,
    где a, b, c, d, f, g принадлежат [0, 1] 
*/

typedef struct gf gf_t;
typedef int gf_inner_t;
typedef struct gf_elem gf_elem_t;

struct gf
{
    int characteristic;     // q
    int power;              // m
    int total_quantity;     // q^m
    gf_inner_t *table;      // таблица элементов поля
    gf_inner_t *rev_table;  // обратная таблица элементов поля
    int forming_polinom;    // полином для построения
    int mask;               // маска для mod операции
};

struct gf_elem
{
    gf_inner_t value;
    gf_t* gf;
};


/*
    Инициализация поля Галуа.
    В текущей реализации q = 2.
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
    Построение таблицы элементов поля.
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
    gf->rev_table = calloc(gf->total_quantity, sizeof(int));

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

    for (int i = 0; i < gf->total_quantity; i++) {
        gf->rev_table[gf->table[i]] = i;
    }
}


/*
    Получение элемента поля по порядковому номеру.
*/

gf_elem_t gf_get(gf_t *gf, int id) {
    gf_elem_t tmp = { gf->table[id], gf };
    return tmp;
}


/*
    Сложение двух элементов поля.
    Реализации для g_elem_t и g_inner_t.
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

gf_inner_t gf_add_inner(gf_inner_t a, gf_inner_t b) {
    return a ^ b;
}


/*
    Умножение двух элементов поля.
    Реализации для g_elem_t и g_inner_t.
*/

gf_elem_t gf_mult(gf_elem_t a, gf_elem_t b) {
    gf_elem_t res = { 0, a.gf };

    if (a.gf != b.gf) {
        res.gf = 0;
        printf("mult error: GF a is not in GF b\n");
    }

    int id1 = a.gf->rev_table[a.value];
    int id2 = a.gf->rev_table[b.value];
    if (id1 == 0 || id2 == 0) res.value = 0;
    else res.value = a.gf->table[((id1+id2-2) % (a.gf->total_quantity-1))+1];
    
    return res;
}

gf_inner_t gf_mult_inner(gf_inner_t a, gf_inner_t b, gf_t *gf) {
    int id1 = gf->rev_table[a];
    int id2 = gf->rev_table[b];
    if (id1 == 0 || id2 == 0) return 0;
    return gf->table[((id1+id2-2) % (gf->total_quantity-1))+1];
}


/*
    Отображение элементов поля через printf.
*/

void gf_print(gf_t *gf) {
    for (int i = 0; i < gf->total_quantity; i++) {
        printf("%d: %d\n", i-1, gf->table[i]);
    }
}


/*
    Очистка памяти, выделенного под поле.
*/

void gf_free(gf_t *gf) {
    free(gf->table);
    free(gf->rev_table);
    free(gf);
}