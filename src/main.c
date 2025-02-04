

#include <stdio.h>
#include "gflib/gflib.h"

int main() {

    int m = 3; // GF(q^m)
    int polinom = 0b1011; // x^3 + x + 1

    gf_t *gf = gf_init(m);
    gf_build(gf, polinom);
    gf_print(gf);
    gf_free(gf);

    return 0;
}