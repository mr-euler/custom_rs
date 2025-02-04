

#include <stdio.h>
#include "gflib/gflib.h"

int main() {

    int m = 3; // GF(q^m)
    int polinom = 0b1011; // x^3 + x + 1

    gf_t *gf = gf_init(m);
    gf_build(gf, polinom);
    // gf_print(gf);

    gf_elem_t a = gf_get(gf, 4);
    gf_elem_t b = gf_get(gf, 5);
    gf_elem_t c = gf_mult(a, b, gf);

    printf("%d * %d = %d\n", a, b, c);

    gf_free(gf);

    return 0;
}