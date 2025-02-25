

#include <stdio.h>
#include "gflib/gflib.h"
#include "gflib/polinom.h"
#include "gflib/gen_polinom.h"

int mypow(int n, int m) {
    int result = 1;
    for (int i = 0; i < m; i++) {
        result *= n;
    }
    return result;
}

int binary_to_decimal(int num) {
    int result = 0;
    int power = 0;
    while (num) {
        // num & 1
        result += mypow(10, power) * (num & 1);
        power++;
        num >>= 1;
    }
    return result;
}

void polynom_print(int result) {
    int first = 1;
    int c = 0;
    while (result) {
        if (result & 1) {
            if (!first) printf(" + ");
            printf("x^%d", c);
            first = 0;
        }
        c++;
        result >>= 1;
    }
    printf("\n");
}

int main() {
    int m; // GF(q^m)
    printf("Введите степень поля Галуа (от 2).\n>>");
    scanf("%d", &m);

    if (m < 3) return 1;

    int from = (1 << m);
    int to = (1 << (m+1));

    for (int i = from; i < to; i++) {
        gf_t *gf = gf_init(m);
        if(!gf_build(gf, i)) {
            printf("%d | ", binary_to_decimal(i));
            polynom_print(i);
            gf_free(gf);
        }
    }
    return 0;
}