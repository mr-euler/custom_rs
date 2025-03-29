
#include <stdio.h>
#include "gflib/gflib.h"
#include "gflib/polinom.h"
#include "gflib/gen_polinom.h"

int count_bits(int num) {
    int count = 0;
    while (num) {
        num >>= 1;
        count++;
    }
    return count;
}

int decimal_to_binary(int num) {
    int decimal = 0;
    int mask = 1;
    while (num) {
        decimal += mask * (num % 10);
        num /= 10;
        mask = mask << 1;
    }
    return decimal;
}

int main() {

    // Параметры кода РС

    // Константы

    int b = 1; // Параметр кодирования

    // Динамические параметры

    // Полином
    // int form_polinom = 0b1011; // x^3 + x + 1
    int form_polinom = 10011;
    printf("Введите образующий полином в 10-тичной форме.\n");
    printf("Пример: x^3+x+1 => 1011\n");
    printf(">>");
    scanf("%d", &form_polinom);
    form_polinom = decimal_to_binary(form_polinom); // Преобразование 1011 => 11
    printf("\n"); // Отступ

    // Расчет параметров поля
    int m = count_bits(form_polinom)-1; // GF(q^m)

    int n = (1 << m)-1; // Количество кодовых символов
    // Подрязумевается, что n+1 == 2^m

    // Количество информационных символов
    int k = 9;
    while (1) {
        printf("Введите количество информационных символов.\nОно должно находиться в пределах от 1 до %d.\n", n-2);
        printf(">>");
        scanf("%d", &k);
        if (1 <= k && k <= n-2) break;
        printf("Неверный ввод.\n");
    }

    // Остальные параметры

    int r = n - k; // Количество исправляющих символов
    int t = r / 2; // Количество возможных для исправления символов 
    int d = n-k+1; // минимальное кодовое расстояние

    // Формирование поля Галуа
    gf_t *gf = gf_init(m);
    if(gf_build(gf, form_polinom)) {
        printf("gf build error: переданный полином не является неприводимым\n");
        return 1;
    }

    // gf_print_bin(gf);

    // 16 2 9 12 2 15 1 2 1 2 1 2 1 2 1
    // gf_elem_t coded[] = { 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 2, 1, 0
        // Владимиров
        // 0,
        // 0,
        // gf_get_by_degree(gf, 1),
        // 0,
        // gf_get_by_degree(gf, 5),
        // 0,
        // 0

        // Студфайл
        // gf_get_by_degree(gf, 5),
        // gf_get_by_degree(gf, 5),
        // gf_get_by_degree(gf, 11),
        // gf_get_by_degree(gf, 12),
        // gf_get_by_degree(gf, 1),
        // gf_get_by_degree(gf, 12),
        // gf_get_by_degree(gf, 13),
        // gf_get_by_degree(gf, 9),
        // gf_get_by_degree(gf, 9),
        // gf_get_by_degree(gf, 14),
        // gf_get_by_degree(gf, 5),
        // gf_get_by_degree(gf, 3),
        // gf_get_by_degree(gf, 2),
        // gf_get_by_degree(gf, 8),
        // gf_get_by_degree(gf, 14)
    // };

    // Чтение данных и проверка на соответствие полю
    gf_elem_t coded[n];
    int read_data = 1;
    while (read_data) {
        read_data = 0;
        printf("Введите полученные кодовые символы (в виде элементов поля) через пробел.\nТаких символов должно быть %d.\n", n);
        printf(">>");
        for (int i = 0; i < n; i++) {
            scanf("%d", &coded[i]);
            if (coded[i] >= gf->total_quantity) read_data = 1;
        }

        if (!read_data) break;
        printf("Получен неверный элемент, повторите попытку.\n");
    }

    polinom_t *encoded = polinom_init(gf, 1);
    polinom_append(encoded, coded, n);
    
    // printf("encoded: ");
    // polinom_print(encoded);
    // printf("encoded size: %d\n", encoded->capacity);

    // Считаем синдромы от b до b+2t-1
    polinom_t *syndroms = polinom_init(gf, b+2*t);
    for (int i = b; i <= b+2*t-1; i++) {
        // printf("%d %d\n", i, gf_get_by_value(gf, polinom_call(encoded, gf_get_by_degree(gf, i))));
        polinom_set(syndroms, i, polinom_call(encoded, gf_get_by_degree(gf, i)));
    }
    // printf("syndrom: ");
    // polinom_print(syndroms);

    // Проверка синдрома на то, содержит ли он только нули
    int syndrom_empty = 1;
    for (int i = 0; i < syndroms->degree; i++) {
        if (syndroms->data[i] != 0) {
            syndrom_empty = 0;
            break;
        }
    }

    if (syndrom_empty) {
        printf("Синдром содержит только нули. Исправление ошибок не требуется\n");
        return 0;
    }

    // Алгоритм БМ

    /*
    // Исходные данные
    polinom_t *sigma = polinom_init(gf, 1);
    polinom_set(sigma, 0, 1);
    polinom_t *B = polinom_init(gf, 1);
    polinom_set(B, 0, 1);
    int L = 0, R = 1, M = 0;
    
    while (R <= d-1) {
        gf_elem_t delta = 0;
        for (int j = 0; j <= L; j++) {
            // Первый аргумент
            gf_elem_t sigma_i;
            if (j >= sigma->capacity) sigma_i = 0;
            else sigma_i = sigma->data[j];
            // Второй аргумент
            gf_elem_t syndroms_i;
            if (R-j-1 >= syndroms->capacity) syndroms_i = 0;
            else syndroms_i = syndroms->data[R-j-1];
            // Перемножение
            delta = gf_add(delta, gf_mult(gf, sigma_i, syndroms_i));
        }

        if (delta == 0) {
            R++;
            continue;
        }

        polinom_t *sigma_new = polinom_copy(sigma); // T(x)
        polinom_t *B_new = polinom_copy(B);
        for (int k = 0; k < B_new->degree; k++) {
            polinom_set(B_new, k, gf_mult(gf, delta, B->data[k]));
        }
        polinom_right_shift(B_new, R-M);
        polinom_add(sigma_new, B_new);
        polinom_free(B_new);

        if (2*L < R) {
            for (int k = 0; k < sigma->degree; k++) {
                polinom_set(sigma, k, gf_div(gf, delta, sigma->data[k]));
            }
            polinom_free(B);
            B = sigma;
            sigma = sigma_new;
            L = R - L;
            M = R;
        } else {
            polinom_free(sigma);
            sigma = sigma_new;
        }
        R++;
    }
    */

    // /*
    // Исходные данные
    // j = b..b+2t-1
    // printf("Метод БМ: начальные условия\n");
    polinom_t *sigma = polinom_init(gf, n+5);
    polinom_set(sigma, 0, 1);
    // printf("sigma: ");
    // polinom_print(sigma);
    polinom_t *p = polinom_init(gf, 2);
    polinom_set(p, 1, 1);
    // printf("p: ");
    // polinom_print(p);
    int l = 0;
    // printf("l = %d\n", l);
    
    for (int i = 1; i < d; i++) {
        // printf("\nИтерация i = %d\n", i);
        // Пункт 1
        gf_elem_t delta = 0;
        for (int j = 0; j <= l; j++) {
            delta = gf_add(delta, gf_mult(gf, sigma->data[j], syndroms->data[i-j-1]));
        }
        // printf("delta = %d", delta);
        // if (delta != 0) printf(" (e^%d)\n", gf_get_by_value(gf, delta));
        // else printf("\n");

        // Пункт 2
        if (delta != 0) {
            // printf("delta != 0\n");
            // Пункт 3
            polinom_t *sigma_new = polinom_copy(sigma);
            polinom_t *p_new = polinom_copy(p);
            for (int k = 0; k < p_new->degree; k++) {
                polinom_set(p_new, k, gf_mult(gf, delta, p->data[k]));
            }
            polinom_add(sigma_new, p_new);
            // printf("1\n");
            polinom_free(p_new);

            // printf("sigma new: ");
            // polinom_print(sigma_new);

            // Пункт 4
            if (2*l < i) {
                // printf("2*l < i\n");
                // Пункт 5
                l = i - l;
                // printf("l = %d\n", l);
                // printf("2\n");
                polinom_free(p);
                p = polinom_copy(sigma);
                for (int k = 0; k < p->degree; k++) {
                    polinom_set(p, k, gf_div(gf, delta, p->data[k]));
                }
                // printf("p: ");
                // polinom_print(p);
            } else {
                // printf("2*l >= i\n");
            }

            // Пункт 6
            // printf("3\n");
            // printf("%d %d %d\n",
            //     sigma->capacity,
            //     sigma->degree,
            //     (sizeof sigma->data) / (sizeof sigma->data[0])
            // );
            polinom_free(sigma);
            sigma = sigma_new;
            // printf("sigma: ");
            // polinom_print(sigma);
        } else {
            // printf("delta == 0\n");
        }

        // Пункт 7
        polinom_right_shift(p, 1);
        // printf("p = xp: ");
        // polinom_print(p);
    }
    // printf("\nКонец\n\n");
    // */

    // gf_elem_t tmp = syndroms->data[syndroms->degree-1];
    // syndroms->data[syndroms->degree-1] = 0;
    // for (int i = syndroms->degree-2; i > -1; i--) {
    //     gf_elem_t tmp2 = syndroms->data[i];
    //     syndroms->data[i] = tmp;
    //     tmp = tmp2;
    // }
    // polinom_calc_degree(syndroms);
    // printf("syndrom new: ");
    // polinom_print(syndroms);

    // Отображение количества ошибок
    if (sigma->degree-1 > t) {
        printf("Количество ошибок в принятой кодовой последовательности превышает допустимый предел исправления (%d) (%d). Декодирование невозможно\n", sigma->degree-1, r);
        return 0;
    }

    printf("Определено количество ошибок: %d\n", sigma->degree-1);

    // w(x) = S(x) * sigma(x) mod x**d
    polinom_t *w = polinom_copy(syndroms);

    // Ультра странный сдвиг, который работает при b = 0
    if (!b) {
        polinom_right_shift(w, 1);
        w->data[0] = 1;
    }

    // printf("w before: ");
    // polinom_print(w);

    polinom_mult(w, sigma);

    // printf("w mult: ");
    // polinom_print(w);

    printf("sigma: ");
    polinom_print(sigma);

    polinom_t *mod = polinom_init(gf, 2*t+1);
    polinom_set(mod, 2*t, 1);
    polinom_mod(w, mod);

    // printf("mod: ");
    // polinom_print(mod);

    // printf("%d\n", mod->capacity);

    polinom_free(mod); // TODO: тут появлялась ошибка

    // TODO: если sigma->degree-1 != t, то выводить информационное сообщение

    // printf("w: ");
    // polinom_print(w);

    // printf("sigma: ");
    // polinom_print(sigma);
    // printf("%d", sigma->degree);
    
    // Определяем количество ошибок
    int error_count = sigma->degree-1;

    // Процедура Ченя (позиции ошибок)
    gf_elem_t solve[error_count]; // содержит элементы Bjl такие, что sigma(Bjl^-1) == 0
    int length = 0;
    for (int i = 0; i < gf->total_quantity; i++) {
        gf_elem_t B = i;
        // printf("B (%d) neg call: %d\n", B, polinom_call(sigma, gf_neg(gf, B)));
        if (polinom_call(sigma, gf_neg(gf, B)) == 0) {
            solve[length] = B;
            length++;
        }
    }

    if (length != error_count) {
        printf("Ошибка при подсчете позиции ошибок (%d)\n", length);
        return 0;
    }

    printf("Позиции ошибок [ ");
    for (int i = 0; i < error_count; i++) {
        printf("x^%d ", gf_get_by_value(gf, solve[i]));
    }
    printf("]\n");

    // Алгоритм форни
    polinom_t *sigma_derivative = polinom_copy(sigma);
    polinom_derivative(sigma_derivative);

    // printf("c'(%d) = %d", solve);
    // printf("sigma derivative: ");
    // polinom_print(sigma_derivative);

    gf_elem_t errors[error_count];
    for (int i = 0; i < error_count; i++) {
        // printf("%d\n", gf_get_by_value(gf, polinom_call(w, gf_neg(gf, solve[i]))));
        errors[i] = gf_mult(
            gf,
            gf_pow(gf, solve[i], 2-b),
            gf_mult(
                gf,
                polinom_call(w, gf_neg(gf, solve[i])),
                gf_neg(gf, polinom_call(sigma_derivative, gf_neg(gf, solve[i])))
            )
        );
    }

    // Отображение ошибок
    printf("Отображение ошибок: [ ");
    for (int i = 0; i < error_count; i++) {
        printf("e^%d ", gf_get_by_value(gf, errors[i]));
    }
    printf("]\n");

    // Проверка
    // printf("encoded (err): ");
    // polinom_print(encoded);

    for (int i = 0; i < error_count; i++) {
        polinom_set(
            encoded,
            gf_get_by_value(gf, solve[i]),
            gf_add(encoded->data[gf_get_by_value(gf, solve[i])], errors[i])
        );
    }

    // printf("encoded: ");
    // polinom_print(encoded);
    // Отображение итоговой кодовой последовательности
    printf("Итоговый полином: [ ");
    for (int i = 0; i < n; i++) {
        printf("%d ", encoded->data[i]);
    }
    printf("]\n");


    printf("Повторная проверка: [ ");
    for (int i = b; i <= b+2*t-1; i++) {
        printf("%d ", polinom_call(encoded, gf_get_by_degree(gf, i)));
    }
    printf("]\n");

    // mod = polinom_copy(encoded);
    // polinom_t *gen_polinom = generating_polinom(gf, b, t);
    // polinom_mod(mod, gen_polinom);
    // printf("Проверка: ");
    // polinom_print(mod);

    gf_free(gf);
    return 0;
}