

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

gf_t* gf_init(int power);
int gf_build(gf_t *gf, int polinom);
gf_elem_t gf_get_by_id(gf_t *gf, int id);
gf_elem_t gf_get_by_degree(gf_t *gf, int id);
gf_elem_t gf_add(gf_elem_t a, gf_elem_t b);
gf_elem_t gf_mult(gf_t *gf, gf_elem_t a, gf_elem_t b);
gf_elem_t gf_pow(gf_t *gf, gf_elem_t elem, int degree);
gf_elem_t gf_div(gf_t *gf, gf_elem_t a, gf_elem_t b);
void gf_print(gf_t *gf);
void gf_free(gf_t *gf);