

typedef struct polinom polinom_t;

struct polinom
{
    gf_t *gf;           // поле Галуа (на котором основан)
    int capacity;       // емкость полинома (количество ячеек памяти для хранения коэффициентов)
    int degree;         // степень полинома
    gf_elem_t *data;    // массив коэффициентов полинома
};

polinom_t* polinom_init(int capacity, gf_t* gf);
void polinom_extencion(polinom_t *polinom, int size);
void polinom_calc_degree(polinom_t *polinom);
void polinom_set(polinom_t *polinom, int index, gf_elem_t elem);
void polinom_mult(polinom_t *polinom1, polinom_t *polinom2);
polinom_t* polinom_copy(polinom_t *polinom1);
void polinom_print(polinom_t *polinom);
void polinom_free(polinom_t *polinom);