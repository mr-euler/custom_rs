

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
    gf_elem_t c = gf_add(a, b);
    gf_elem_t d = gf_mult(a, b);

    printf("%d + %d = %d\n", a.value, b.value, c.value);
    printf("%d * %d = %d\n", a.value, b.value, d.value);

    gf_free(gf);

    return 0;
}