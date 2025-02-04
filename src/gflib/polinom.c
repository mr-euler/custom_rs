
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gflib.h"

typedef struct polinom polinom_t;

struct polinom
{
    gf_t *gf;
    int capacity;
    int degree;
    gf_inner_t* data;
};

polinom_t* polinom_init(int capacity, gf_t* gf) {
    polinom_t *polinom = malloc(sizeof(polinom_t));
    polinom->gf = gf;
    polinom->capacity = capacity;
    polinom->degree = 0;
    polinom->data = calloc(capacity, sizeof(gf_inner_t));
    for (int i = 0; i < capacity; i++) polinom->data[i] = 0;
    return polinom;
}

void polinom_extencion(polinom_t *polinom) {
    gf_inner_t* new_data = calloc(polinom->capacity * 2, sizeof(gf_inner_t));
    memcpy(new_data, polinom->data, polinom->capacity * sizeof(gf_inner_t));
    free(polinom->data);
    polinom->data = new_data;
    polinom->capacity *= 2;
}

void polinom_calc_degree(polinom_t *polinom) {
    polinom->degree = 0;
    for (int i = 0; i < polinom->capacity; i++) {
        if (polinom->data[i] != 0) polinom->degree = i+1;
    }
}

void polinom_set(polinom_t *polinom, int index, gf_elem_t elem) {
    
    if (polinom->gf != elem.gf) {
        printf("set error: polinom is not in GF\n");
        return;
    }

    if (index >= polinom->capacity) {
        printf("set error: index is out of capaciry\n");
    }

    polinom->data[index] = elem.value;
    if (polinom->degree <= index) polinom->degree = index+1;
    if (elem.value == 0) polinom_calc_degree(polinom);
}

void polinom_set_inner(polinom_t *polinom, int index, gf_inner_t elem) {

    if (index >= polinom->capacity) {
        printf("set error: index is out of capaciry\n");
    }

    polinom->data[index] = elem;
    if (polinom->degree <= index) polinom->degree = index+1;
    if (elem == 0) polinom_calc_degree(polinom);
}


polinom_t* polinom_mult(polinom_t *polinom1, polinom_t *polinom2) {
    polinom_t *polinom = polinom_init(polinom1->capacity + polinom2->capacity, polinom1->gf);

    if (polinom1->gf != polinom2->gf) {
        printf("polinom mult error: polinoms are not in same GF\n");
        return polinom;
    }

    for (int i = 0; i < polinom1->degree; i++) {
        for (int j = 0; j < polinom2->degree; j++) {
            gf_inner_t mult = gf_mult_inner(polinom1->data[i], polinom2->data[j], polinom1->gf);
            polinom_set_inner(polinom, i+j, gf_add_inner(polinom->data[i+j], mult));
        }
    }

    return polinom;

}



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

void polinom_free(polinom_t *polinom) {
    free(polinom->data);
    free(polinom);
}