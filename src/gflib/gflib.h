

typedef struct gf gf_t;
typedef int gf_inner_t;
typedef struct gf_elem gf_elem_t;

struct gf
{
    int characteristic;     // q
    int power;              // m
    int total_quantity;     // q^m
    gf_inner_t *table;      // таблица элементов поля
    gf_inner_t *rev_table;  // обратная таблица элементов поля
    int forming_polinom;    // полином для построения
    int mask;               // маска для mod операции
};

struct gf_elem
{
    gf_inner_t value;
    gf_t* gf;
};

gf_t* gf_init(int power);
int gf_build(gf_t *gf, int polinom);
gf_elem_t gf_get(gf_t *gf, int id);
gf_elem_t gf_add(gf_elem_t a, gf_elem_t b);
gf_inner_t gf_add_inner(gf_inner_t a, gf_inner_t b);
gf_elem_t gf_mult(gf_elem_t a, gf_elem_t b);
gf_inner_t gf_mult_inner(gf_inner_t a, gf_inner_t b, gf_t *gf);
void gf_print(gf_t *gf);
void gf_free(gf_t *gf);