
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
typedef int gf_elem_t;

struct gf
{
    int characteristic;     // q
    int power;              // m
    int total_quantity;     // q^m
    gf_elem_t *table;      // таблица элементов поля
    gf_elem_t *rev_table;  // обратная таблица элементов поля
    int forming_polinom;    // полином для построения
    int mask;               // маска для mod операции
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

int gf_build(gf_t *gf, int polinom) {
    gf->table = calloc(gf->total_quantity, sizeof(int));
    gf->rev_table = calloc(gf->total_quantity, sizeof(int));
    
    // Специальный массив для проверки приводимости полинома
    gf_elem_t *counter = calloc(gf->total_quantity, sizeof(int));
    for (int i = 0; i < gf->total_quantity; i++) {
        counter[i] = 0;
    }


    gf->forming_polinom = polinom;

    gf->mask = 1 << gf->power;
    gf->table[0] = 0;
    gf->table[1] = 1;

    counter[0] = 1;
    counter[1] = 1;

    for (int i = 2; i < gf->total_quantity; i++) {
        gf->table[i] = gf->table[i-1] << 1;
        if(gf->table[i] & gf->mask) {
            gf->table[i] ^= polinom;
        }
        
        counter[gf->table[i]]++;
        if (counter[gf->table[i]] > 1) {
            free(gf->table);
            free(gf->rev_table);
            free(counter);
            return 1;
        }
    }

    for (int i = 0; i < gf->total_quantity; i++) {
        gf->rev_table[gf->table[i]] = i;
    }

    return 0;
}


/*
    Получение элемента поля по порядковому номеру примитивного элемента.
*/

gf_elem_t gf_get_by_id(gf_t *gf, int id) {
    if (id < 0) {
        printf("gf get error: invalid index {%d}", id);
        return 0;
    }
    return gf->table[id % gf->total_quantity];
}

/*
    Получение элемента поля по степени примитивного элемента.
*/

gf_elem_t gf_get_by_degree(gf_t *gf, int degree) {
    if (degree < 0) {
        printf("gf get degree error: invalid index {%d}", degree);
        return 1;
    }
    return gf->table[(degree % (gf->total_quantity-1))+1];
}


/*
    Сложение двух элементов поля.
    Реализации для g_elem_t и g_inner_t.
*/

gf_elem_t gf_add(gf_elem_t a, gf_elem_t b) {
    return a ^ b;
}


/*
    Умножение двух элементов поля.
    Реализации для g_elem_t и g_inner_t.
*/

gf_elem_t gf_mult(gf_t *gf, gf_elem_t a, gf_elem_t b) {
    int id1 = gf->rev_table[a];
    int id2 = gf->rev_table[b];
    if (id1 == 0 || id2 == 0) return 0;
    return gf->table[((id1+id2-2) % (gf->total_quantity-1))+1];
}


/*
    Возведение элемента поля в степень.
*/

gf_elem_t gf_pow(gf_t *gf, gf_elem_t elem, int degree) {
    if (elem == 0) return 0;
    return gf->table[(( (gf->rev_table[elem]-1) * degree ) % ( gf->total_quantity-1))+1];
}


/*
    Нахождение такого элемента поля,
    который нужно умножить на первый
    аргумент, чтобы получить второй.
    (деление)
*/

gf_elem_t gf_div(gf_t *gf, gf_elem_t a, gf_elem_t b) {
    if (a == 0 && b == 0) return 0;
    if (a == 0 || b == 0) {
        printf("gf div error: invalid args\n");
        return 0;
    }
    int elem_a_degree = gf->rev_table[a]-1;
    int elem_b_degree = gf->rev_table[b]-1;
    if (elem_b_degree >= elem_a_degree) return gf->table[elem_a_degree + elem_b_degree + 1];
    else return gf->table[(gf->total_quantity-1) - elem_a_degree + elem_b_degree + 1];
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