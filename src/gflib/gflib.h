

struct gf
{
    int characteristic;     // q
    int power;              // m
    int total_quantity;     // q^m
    int *table;             // таблица элементов поля
    int polinom;            // полином для построения
};

typedef struct gf gf_t;
typedef int gf_elem_t;

gf_t* gf_init(int power);
void gf_build(gf_t *gf, int polinom);
void gf_print(gf_t *gf);
void gf_free(gf_t *gf);