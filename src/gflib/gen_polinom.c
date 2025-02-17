
#include "gflib.h"
#include "polinom.h"

/*
    Функция для формирования порождающего полинома
*/

polinom_t* generating_polinom(gf_t *gf, int b, int t) {
    polinom_t *polinom_res = polinom_init(gf, 2);
    polinom_set(polinom_res, 0, gf_get_by_id(gf, 1+b));
    polinom_set(polinom_res, 1, gf_get_by_id(gf, 1));
    for (int i = b+1; i < b+2*t; i++) {
        polinom_t *tmp = polinom_init(gf, 2);
        polinom_set(tmp, 0, gf_get_by_id(gf, 1+i));
        polinom_set(tmp, 1, gf_get_by_id(gf, 1));

        polinom_mult(polinom_res, tmp);
        polinom_free(tmp);
    }
    return polinom_res;
}