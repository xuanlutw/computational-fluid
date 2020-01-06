# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include <math.h>

# define N   20
# define DIM 3
# define DX  1.0
# define EPS 0.0000000001
# define GAUSS_ITER_MAX  3000
# define ITER_MAX 5000
# define TEST_SIZE 20
# define DT        0.01

typedef struct {
    int grid_size;
    double* val;
} Field_s;

typedef struct {
    Field_s* comp[DIM];
} Field_v;

typedef struct {
    int num;
    int cent;
    int r;
    int*    coord[DIM];
    double* da[DIM];
} Ball_pts;

typedef struct {
    double val[DIM];
} Vector;

Field_s* init_s(int grid_size) {
    Field_s* s   = (Field_s*) malloc(sizeof(Field_s));
    s->grid_size = grid_size;
    s->val       = malloc(sizeof(double) * grid_size * grid_size * grid_size);
    return s;
}
# define Val_s(P, i, j, k) (*((P)->val + (((i) * (P)->grid_size + (j)) * (P)->grid_size) + (k)))

Field_v* init_v(int grid_size) {
    Field_v* v   = (Field_v*) malloc(sizeof(Field_v));
    for (int i = 0; i < DIM; ++i)
        v->comp[i] = init_s(grid_size);
    return v;
}
# define Val_v(V, s, i, j, k) (Val_s((V)->comp[(s)], i, j, k))

void Lap_s (Field_s* src, Field_s* dest) {
    // dest = nabla^2 src
    assert(src->grid_size == dest->grid_size);
    // skip out boundary
    for (int i = 1; i < src->grid_size - 1; ++i) {
        for (int j = 1; j < src->grid_size - 1; ++j) {
            for (int k = 1; k < src->grid_size - 1; ++k) {
                Val_s(dest, i, j, k) = -6.0 * Val_s(src, i, j, k);
                Val_s(dest, i, j, k) += Val_s(src, i + 1, j, k);
                Val_s(dest, i, j, k) += Val_s(src, i - 1, j, k);
                Val_s(dest, i, j, k) += Val_s(src, i, j + 1, k);
                Val_s(dest, i, j, k) += Val_s(src, i, j - 1, k);
                Val_s(dest, i, j, k) += Val_s(src, i, j, k + 1);
                Val_s(dest, i, j, k) += Val_s(src, i, j, k - 1);
                Val_s(dest, i, j, k) /= DX * DX;
            }
        }
    }
    return;
}

void Lap_v (Field_v* src, Field_v* dest) {
    for (int i = 0; i < DIM; ++i)
        Lap_s(src->comp[i], dest->comp[i]);
    return;
}

void test_Lap_s() {
    Field_s* src  = init_s(TEST_SIZE);
    Field_s* dest = init_s(TEST_SIZE);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k)
                Val_s(src, i, j, k) = i * i + j * j + k * k;
    Lap_s(src, dest);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k)
                printf("%d %d %d %lf\n", i, j, k, Val_s(dest, i, j, k));
    return;
}

void Div(Field_v* src, Field_s* dest){
    // dest = nabla dot src
    assert(src->comp[0]->grid_size == dest->grid_size);
    // skip out boundary
    for (int i = 1; i < dest->grid_size - 1; ++i) {
        for (int j = 1; j < dest->grid_size - 1; ++j) {
            for (int k = 1; k < dest->grid_size - 1; ++k) {
                Val_s(dest, i, j, k) = 0;
                Val_s(dest, i, j, k) += Val_v(src, 0, i + 1, j, k);
                Val_s(dest, i, j, k) -= Val_v(src, 0, i - 1, j, k);
                Val_s(dest, i, j, k) += Val_v(src, 1, i, j + 1, k);
                Val_s(dest, i, j, k) -= Val_v(src, 1, i, j - 1, k);
                Val_s(dest, i, j, k) += Val_v(src, 2, i, j, k + 1);
                Val_s(dest, i, j, k) -= Val_v(src, 2, i, j, k - 1);
                Val_s(dest, i, j, k) /= 2 * DX;
            }
        }
    }
    return;
}

void test_Div() {
    Field_v* src  = init_v(TEST_SIZE);
    Field_s* dest = init_s(TEST_SIZE);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k) {
                Val_v(src, 0, i, j, k) = i;
                Val_v(src, 1, i, j, k) = j;
                Val_v(src, 2, i, j, k) = k;
            }
    Div(src, dest);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k)
                printf("%d %d %d %lf\n", i, j, k, Val_s(dest, i, j, k));
    return;
}

void Sub_D(Field_v* src1, Field_v* src2, Field_v* dest) {
    // dest = (src1 dot nabla) src2
    assert(src1->comp[0]->grid_size == dest->comp[0]->grid_size);
    assert(src2->comp[0]->grid_size == dest->comp[0]->grid_size);
    // skip out boundary
    for (int i = 1; i < dest->comp[0]->grid_size - 1; ++i) {
        for (int j = 1; j < dest->comp[0]->grid_size - 1; ++j) {
            for (int k = 1; k < dest->comp[0]->grid_size - 1; ++k) {
                for (int s = 0; s < DIM; ++s) {
                    Val_v(dest, s, i, j, k) = 0;
                    Val_v(dest, s, i, j, k) += (Val_v(src2, s, i + 1, j, k) - Val_v(src2, s, i - 1, j, k)) * Val_v(src1, 0, i, j, k) / (2 * DX);
                    Val_v(dest, s, i, j, k) += (Val_v(src2, s, i, j + 1, k) - Val_v(src2, s, i, j - 1, k)) * Val_v(src1, 1, i, j, k) / (2 * DX);
                    Val_v(dest, s, i, j, k) += (Val_v(src2, s, i, j, k + 1) - Val_v(src2, s, i, j, k - 1)) * Val_v(src1, 2, i, j, k) / (2 * DX);
                }
            }
        }
    }
    return;
}

void test_Sub_D() {
    Field_v* src  = init_v(TEST_SIZE);
    Field_v* dest = init_v(TEST_SIZE);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k) {
                Val_v(src, 0, i, j, k) = i;
                Val_v(src, 1, i, j, k) = 0;
                Val_v(src, 2, i, j, k) = 0;
            }
    Sub_D(src, src, dest);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k)
                printf("%d %d %d %lf %lf %lf\n", i, j, k, Val_v(dest, 0, i, j, k), Val_v(dest, 1, i, j, k), Val_v(dest, 2, i, j, k));
    return;
}

void Gauss(Field_s* src, Field_s* dest) {
    // dest = nabla^2 src
    assert(src->grid_size == dest->grid_size);
    for (int iter = 0; iter < GAUSS_ITER_MAX; ++iter) {
        double diff = 0.0;
        //printf("%d\n", iter);
        for (int i = 1; i < dest->grid_size - 1; ++i) {
            for (int j = 1; j < dest->grid_size - 1; ++j) {
                for (int k = 1; k < dest->grid_size - 1; ++k) {
                    double tmp = 0;
                    tmp += Val_s(dest, i + 1, j, k);
                    tmp += Val_s(dest, i - 1, j, k);
                    tmp += Val_s(dest, i, j + 1, k);
                    tmp += Val_s(dest, i, j - 1, k);
                    tmp += Val_s(dest, i, j, k + 1);
                    tmp += Val_s(dest, i, j, k - 1);
                    tmp -= Val_s(src, i, j, k) * DX * DX;
                    tmp /= 6.0;
                    diff += (Val_s(dest, i, j, k) - tmp) * (Val_s(dest, i, j, k) - tmp);
                    Val_s(dest, i, j, k) = tmp;
                }
            }
        }
        if (diff < EPS)
            return;
    }
}

void test_Gauss() {
    Field_s* src  = init_s(TEST_SIZE);
    Field_s* dest = init_s(TEST_SIZE);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k) {
                if ((i == TEST_SIZE / 2) && (i == j) && (j == k))
                    Val_s(src, i, j, k) = 1;
                else
                    Val_s(src, i, j, k) = 0;
            }
    Gauss(src, dest);
    for (int i = 0; i < TEST_SIZE; ++i)
        printf("%d %lf\n", i, Val_s(dest, i, TEST_SIZE / 2, TEST_SIZE / 2));
    return;
}

void Grad(Field_s* src, Field_v* dest) {
    // dest = nabla src
    assert(src->grid_size == dest->comp[0]->grid_size);
    // skip out boundary
    for (int i = 1; i < src->grid_size - 1; ++i) {
        for (int j = 1; j < src->grid_size - 1; ++j) {
            for (int k = 1; k < src->grid_size - 1; ++k) {
                Val_v(dest, 0, i, j, k) = (Val_s(src, i + 1, j, k) - Val_s(src, i - 1, j, k)) / (2 * DX);
                Val_v(dest, 1, i, j, k) = (Val_s(src, i, j + 1, k) - Val_s(src, i, j - 1, k)) / (2 * DX);
                Val_v(dest, 2, i, j, k) = (Val_s(src, i, j, k + 1) - Val_s(src, i, j, k - 1)) / (2 * DX);
            }
        }
    }
    return;
}

void test_Grad() {
    Field_s* src  = init_s(TEST_SIZE);
    Field_v* dest = init_v(TEST_SIZE);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k) 
                Val_s(src, i, j, k) = i;
    Grad(src, dest);
    for (int i = 0; i < TEST_SIZE; ++i)
        for (int j = 0; j < TEST_SIZE; ++j)
            for (int k = 0; k < TEST_SIZE; ++k)
                printf("%d %d %d %lf %lf %lf\n", i,j, k, Val_v(dest, 0, i, j, k), Val_v(dest, 1, i, j, k), Val_v(dest, 2, i, j, k));
    return;
}

void Arith_s(Field_s* src1, Field_s* src2, Field_s* dest, double coef1, double coef2) {
    // dest = coef1 * src1 + coef2 * src2
    assert(src1->grid_size == dest->grid_size);
    assert(src2->grid_size == dest->grid_size);

    for (int i = 0; i < src1->grid_size; ++i) {
        for (int j = 0; j < src1->grid_size; ++j) {
            for (int k = 0; k < src1->grid_size; ++k) {
                Val_s(dest, i, j, k) = coef1 * Val_s(src1, i, j, k) + coef2 * Val_s(src2, i, j, k);
            }
        }
    }
    return;
}

void Arith_v (Field_v* src1, Field_v* src2, Field_v* dest, double coef1, double coef2) {
    // dest = coef1 * src1 + coef2 * src2
    for (int i = 0; i < DIM; ++i)
        Arith_s(src1->comp[i], src2->comp[i], dest->comp[i], coef1, coef2);
    return;
}

void inner_bound_rect(Field_v* src, int l) {
    // Rect
    int center = src->comp[0]->grid_size / 2;
    for (int i = center - l / 2; i < center + l / 2; ++i)
        for (int j = center - l / 2; j < center + l / 2; ++j)
            for (int k = center - l / 2; k < center + l / 2; ++k)
                for (int s = 0; s < DIM; ++s)
                    Val_v(src, s, i, j, k) = 0;
    return;
}

void inner_bound_ball(Field_v* src, int r) {
    int g_size = src->comp[0]->grid_size;
    int cent = src->comp[0]->grid_size / 2;
    for (int i = 0; i < g_size; ++i)
        for (int j = 0; j < g_size; ++j)
            for (int k = 0; k < g_size; ++k)
                if ((i - cent) * (i - cent) + (j - cent) * (j - cent) + (k - cent) * (k - cent) <= r * r)
                    for (int s = 0; s < DIM; ++s)
                        Val_v(src, s, i, j, k) = 0;
    return;
}

Ball_pts* init_ball_pts(int r, int size) {
    Ball_pts* b = malloc(sizeof(Ball_pts));
    int pts_max = (int)(400 * M_PI * r * r);
    int cent = size / 2;
    b->r    = r;
    b->num  = 0;
    b->cent = cent;
    for (int i = 0; i < DIM; ++i) {
        b->coord[i] = malloc(sizeof(int) * pts_max);
        b->da[i] = malloc(sizeof(double) * pts_max);
    }
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            for (int k = 0; k < size; ++k) {
                int dist_2 = (i - cent) * (i - cent) + (j - cent) * (j - cent) + (k - cent) * (k - cent);
                double dist = sqrt(dist_2);
                if (dist_2 >= r * r && dist_2 < (r + 1) * (r + 1)) {
                    b->coord[0][b->num] = i;
                    b->coord[1][b->num] = j;
                    b->coord[2][b->num] = k;
                    b->da[0][b->num] = (double)(i - cent) / r;
                    b->da[1][b->num] = (double)(j - cent) / r;
                    b->da[2][b->num] = (double)(k - cent) / r;
                    ++(b->num);
                }
            }
    return b;
}

Vector Comp_force_ball (Field_v* v, Field_s* p, double mu, Ball_pts* pts) {
    Vector force;
    for (int i = 0; i < DIM; ++i) 
        force.val[i] = 0;
    for (int n = 0; n < pts->num; ++n) {
        int x = pts->coord[0][n];
        int y = pts->coord[1][n];
        int z = pts->coord[2][n];
        double dax = pts->da[0][n];
        double day = pts->da[1][n];
        double daz = pts->da[2][n];

        force.val[0] += -Val_s(p, x, y, z) * dax;
        force.val[1] += -Val_s(p, x, y, z) * day;
        force.val[2] += -Val_s(p, x, y, z) * daz;

        force.val[0] += (Val_v(v, 0, x + 1, y, z) - Val_v(v, 0, x - 1, y, z)) / DX * dax;
        force.val[1] += (Val_v(v, 1, x, y + 1, z) - Val_v(v, 1, x, y - 1, z)) / DX * day;
        force.val[2] += (Val_v(v, 2, x, y, z + 1) - Val_v(v, 2, x, y, z - 1)) / DX * daz;

        force.val[0] += (Val_v(v, 0, x, y + 1, z) - Val_v(v, 0, x, y - 1, z)) / 2.0 / DX * day;
        force.val[0] += (Val_v(v, 0, x, y, z + 1) - Val_v(v, 0, x, y, z - 1)) / 2.0 / DX * daz;
        force.val[1] += (Val_v(v, 1, x + 1, y, z) - Val_v(v, 1, x - 1, y, z)) / 2.0 / DX * dax;
        force.val[1] += (Val_v(v, 1, x, y, z + 1) - Val_v(v, 1, x, y, z - 1)) / 2.0 / DX * daz;
        force.val[2] += (Val_v(v, 2, x + 1, y, z) - Val_v(v, 2, x - 1, y, z)) / 2.0 / DX * dax;
        force.val[2] += (Val_v(v, 2, x, y + 1, z) - Val_v(v, 2, x, y - 1, z)) / 2.0 / DX * day;

        force.val[0] += (Val_v(v, 1, x + 1, y, z) - Val_v(v, 1, x - 1, y, z)) / 2.0 / DX * day;
        force.val[0] += (Val_v(v, 2, x + 1, y, z) - Val_v(v, 2, x - 1, y, z)) / 2.0 / DX * daz;
        force.val[1] += (Val_v(v, 0, x, y + 1, z) - Val_v(v, 0, x, y - 1, z)) / 2.0 / DX * dax;
        force.val[1] += (Val_v(v, 2, x, y + 1, z) - Val_v(v, 2, x, y - 1, z)) / 2.0 / DX * daz;
        force.val[2] += (Val_v(v, 0, x, y, z + 1) - Val_v(v, 0, x, y, z - 1)) / 2.0 / DX * dax;
        force.val[2] += (Val_v(v, 1, x, y, z + 1) - Val_v(v, 1, x, y, z - 1)) / 2.0 / DX * day;
    }
    for (int i = 0; i < DIM; ++i) 
        force.val[i] *= (4.0 * M_PI * pts->r * pts->r / pts->num);
    return force;
}

int main () {
    //test_Lap_s();
    //test_Div();
    //test_Sub_D();
    //test_Gauss();
    //test_Grad();
    
    Field_v* v[2];
    v[0] = init_v(N);
    v[1] = init_v(N);

    // IC
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                for (int s = 0; s < DIM; ++s) {
                    if ((i == N - 1 || j == N - 1 || k == N - 1 || i == 0 || j == 0 || k == 0) && s == 0) {
                        Val_v(v[0], s, i, j, k) = 1;
                        Val_v(v[1], s, i, j, k) = 1;
                    }
                    else {
                        Val_v(v[0], s, i, j, k) = 0;
                        Val_v(v[1], s, i, j, k) = 0;
                    }
                }
            }
        }
    }
    double rho = 1.0;
    double mu  = 1.;

    Field_v* res1      = init_v(N);
    Field_v* res2      = init_v(N);
    Field_v* pre_p_src = init_v(N);
    Field_s* p_src     = init_s(N);
    Field_s* p         = init_s(N);
    Field_v* dp        = init_v(N);
    Field_v* dv        = init_v(N);
    Ball_pts* pts      = init_ball_pts(4, N);

    for (int iter = 0; iter < ITER_MAX; ++iter) {
        int idx_now = iter % 2;
        int idx_nxt = (iter + 1) % 2;
        
        Lap_v(v[idx_now], res1);
        Sub_D(v[idx_now], v[idx_now], res2);
        Arith_v(res1, res2, pre_p_src, mu, -rho);
        Div(pre_p_src, p_src);
        Gauss(p_src, p);
        Grad(p, dp);
        Arith_v(dp, pre_p_src, dv, -1.0 / rho, 1.0 / rho);
        Arith_v(v[idx_now], dv, v[idx_nxt], 1.0, DT);
        //inner_bound_rect(v[idx_nxt], 6);
        inner_bound_ball(v[idx_nxt], 4);
        Vector force = Comp_force_ball(v[idx_nxt], p, mu, pts);
        
        // print
#define erase_display_all() printf("\033[2J");fflush(stdout)
        erase_display_all();
        printf("%d\n", iter);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                printf("%.3lf ", Val_v(v[idx_nxt], 0, i, j, N / 2));
                /*
                if (i >= 7 && i < 13 && j >= 7 && j < 13)
                    //printf("\033[31m%.3lf ", Val_s(p, i, N / 2, j));
                    printf("\033[31m%.3lf ", Val_v(v[idx_nxt], 0, i, N / 2, j));
                else
                    //printf("\033[37m%.3lf ", Val_s(p, i, N / 2, j));
                    printf("\033[37m%.3lf ", Val_v(v[idx_nxt], 0, i, N / 2, j));
                    */
                /*
                if ((i - N / 2) * (i - N / 2) + (j - N / 2) * (j - N / 2) <= 25)
                    printf("\033[31m%.3lf ", Val_v(v[0], 1, i, j, N / 2));
                else
                    printf("\033[37m%.3lf ", Val_v(v[0], 1, i, j, N / 2));
                    */
            }
            printf("\n");
        }
        printf("\n%lf %lf %lf\n", force.val[0], force.val[1], force.val[2]);
    }

    return 0;
}
