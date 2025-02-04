

typedef struct gf gf_t;
typedef int gf_elem_t;

struct gf
{
    int characteristic;     // q
    int power;              // m
    int total_quantity;     // q^m
    gf_elem_t *table;       // таблица элементов поля
    int polinom;            // полином для построения
    int mask;               // маска для mod операции
};

gf_t* gf_init(int power);
void gf_build(gf_t *gf, int polinom);
gf_elem_t gf_get(gf_t *gf, int id);
gf_elem_t gf_add(gf_elem_t a, gf_elem_t b);
gf_elem_t gf_mult(gf_elem_t a, gf_elem_t b, gf_t *gf);
void gf_print(gf_t *gf);
void gf_free(gf_t *gf);