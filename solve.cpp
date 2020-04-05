#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/sysinfo.h>

#include "mpi.h"
#include "matrix.h"
#include "solve.h"

static double eps = 1e-14;

void get_block(double *A, double *B, int n, int m, int l, int i, int j, int real)
{
    int p, q = 0;
    int P = real * (i - 1) * n + real * (j - 1), pp = P;
    int hm = P + (m - 1) * n;

    while (pp <= hm)
    {
        for (p = pp; p < pp + l; p++)
        {
            B[q] = A[p];
            q = q + 1;
        }
        q += real - l;
        pp += n;
    }
}

inline void put_block(double *A, double *B, int n, int m, int l, int i, int j, int real)
{
    int p, q = 0;
    int P = real * (i - 1) * n + real * (j - 1), pp = P;
    int hm = P + (m - 1) * n;

    while (pp <= hm)
    {
        for (p = pp; p < pp + l; p++)
        {
            A[p] = B[q];
            q = q + 1;
        }
        q += real - l;
        pp += n;
    }
}

void print(double *a, int m, int l, int real)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < l; j++)
            printf("%.5f ",a[i*real + j]);
        printf("\n");
    }
    printf("\n");
}

/* блочное умножение блоков */
inline void block_mult_block(double *a, double *b, double *c, int m, int n, int l, int real) // (m x n) * (n x l)
{
    double s00 = 0, s01 = 0, s10 = 0, s02 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;
    double *pa, *pb, *pc;
    double *pa1, *pa2;
    double *pc1, *pc2;
    double b0, b1, b2;
    double a0, a1, a2;
    int i, j, k;

    for (i = 0; i < m; i += 3)
    {
        pa = a + i * real;
        pa1 = a + (i + 1) * real;
        pa2 = a + (i + 2) * real;

        pc = c + i * real;
        pc1 = c + (i + 1) * real;
        pc2 = c + (i + 2) * real;

        for (j = 0; j < l; j += 3)
        {
            s00 = 0; s01 = 0; s02 = 0;
            s10 = 0; s11 = 0; s12 = 0;
            s20 = 0; s21 = 0; s22 = 0.;
            for (k = 0; k < n; ++k)
            {
                pb = b + k * real + j;

                b0 = pb[0];
                b1 = pb[1];
                b2 = pb[2];

                a0 = pa[k];
                a1 = pa1[k];
                a2 = pa2[k];

                s00 += a0 * b0;
                s01 += a1 * b0;
                s02 += a2 * b0;

                s10 += a0 * b1;
                s11 += a1 * b1;
                s12 += a2 * b1;

                s20 += a0 * b2;
                s21 += a1 * b2;
                s22 += a2 * b2;
            }
            pc[j] = s00;
            pc1[j] = s01;
            pc2[j] = s02;

            pc[j+1] = s10;
            pc1[j+1] = s11;
            pc2[j+1] = s12;

            pc[j+2] = s20;
            pc1[j+2] = s21;
            pc2[j+2] = s22;
        }
    }
}

/* блок размера(m x n) - блок размера(m x n) */
inline void block_minus_block(double *a, double *b, double *c, int m, int n, int real)
{
    int i, j;
    double *pc;
        for (i = 0; i < m; ++i)
        {
            pc = c + i * real;
            for (j = 0; j < n; ++j)
                pc[j] = a[i * real + j] - b[i * real + j];
        }
}

/* блок размера(m x n) + блок размера(m x n) */
inline void block_plus_block(double *a, double *b, double *c, int m, int n, int real)
{
    int i, j;
    double *pc;
        for (i = 0; i < m; ++i)
        {
            pc = c + i * real;
            for (j = 0; j < n; ++j)
                pc[j] = a[i * real + j] + b[i * real + j];
        }
}

/* положить в блок единичную матрицу */
inline void put_block_E(double *a, int m, int l)// A(m x l) = E
{
    int i, j;

    for (i = 0; i < m; ++i)
        for (j = 0; j < l; ++j)
        {
            if (j != i) a[i * l + j] = 0;
            else a[i * l + j] = 1;
        }
}

/* поделить элементы i-той строки на (i,i)-тый элемент */
inline int division1(double *a, double *e, int i, int n, double z, int m, double z1)
{
    int j;
    double f = a[i * m + i];

    if ((fabs(f) < eps * z) && (n > 1)) return -1;
    else if ((n == 1) && (fabs(f) < eps * z1)) return -1;
    else
    {
        a[i * m + i] = 1;
        /* для матрицы А */
        for (j = i + 1; j < n; ++j)
          a[i * m + j] /= f;
          /* для присоединенной матрицы */
        for (j = 0; j <= i; ++j)
          e[i * m + j] /= f;
    }
    return 0;
}

/* номер первого столбца, с которым начинает работу q-тый поток в i-той строке */
int f(int i, int q, int p)
{
    int k = i / p;
    if ((i % p) > q) return k * p + q + p;
    else return k * p + q;
}

/* найти обратный к (i,i)-тому блоку, умножить блоки i-той строки на обратный к (i,i)-тому блоку */
int division(double *a, double *e, double *x, double *y, double *z, int i, int n, int m, int l, int k, int q, int p, double *h, double z1)
{
    int j, w, ww, g, s;

    g = (l == 0 ? k : k + 1);

    if ((i == k + 1) && (l != 0)) w = l;
    else w = m;

    get_block(h + (q - 1) * n * m, x + (q - 1) * m * m, m, w, w, i, 1, m);

    put_block_E(y + (q - 1) * m * m, m, m);

    /* находим обратную */
    s = solve1(x + (q - 1) * m * m, y + (q - 1) * m * m, w, m, z1);

    if (s < 0) return -1;
    else
    {
        int ff = f(i, q, p);
        /* умножаем на обратную i-тую строку матрицы A*/
        for (j = ff; j <= g; j += p)
        {
            if (j > i)
            {
                if ((j == k + 1) && (l != 0)) ww = l;
                else ww = m;
                get_block(a, x + (q - 1) * m * m, n, w, ww, i, j, m);
                block_mult_block(y + (q - 1) * m * m, x + (q - 1) * m * m, z + (q - 1) * m * m, w, w, ww, m);
                put_block(a, z + (q - 1) * m * m, n, w, ww, i, j, m);
            }
        }
        /* умножаем на обратную i-тую строку присоединенной матрицы*/
        for (j = q; j <= i; j += p)
        {
            if ((j == g) && (l != 0)) ww = l;
            else ww = m;
            get_block(e, x + (q - 1) * m * m, n, w, ww, i, j, m);
            block_mult_block(y + (q - 1) * m * m, x + (q - 1) * m * m, z + (q - 1) * m * m, w, w, ww, m);
            put_block(e, z + (q - 1) * m * m, n, w, ww, i, j, m);
        }
    }
    return 0;
}

/* цикл вычитания для простого варианта */
inline void subtraction1(double *a, double *e, int i, int n, int m)
{
    int str, j;
    double f;

    for (str = i + 1; str < n; ++str)
    {
        f = a[str * m + i];
        a[str * m + i] = 0;
        /* для матрицы А */
        for (j = i + 1; j < n; ++j)
            a[str * m + j] -= f * a[i * m + j];
        /* для присоединенной матрицы */
        for (j = 0; j < str; ++j)
            e[str * m + j] -= f * e[i * m + j];
    }
}

/* цикл вычитания для MPI варианта */
void subtraction(double *a, double *b, double *x, double *y, double *z, double *h, double *buf, int n, int m, int l, int k, int max_cols, int i, int local_ff)
{
    int q, p, size, size1;
    int s, j, global_j;
    double *aa;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    for (j = local_ff; j <= max_cols; j++)
    {
        global_j = (j - 1) * p + q + 1;

        if (j == max_cols && l > 0 && global_j == k) size1 = l;
        else size1 = m;

        aa = a + (j - 1) * n * m;

        for (s = i + 1; s <= k; ++s)
        {
            if (s != k) size = m;
            else if (l > 0) size = l;
            else size = m;

            // i-тый столбец уже лежит в buf
            get_block(buf, h, m, size, m, s, 1, m);
            get_block(aa, y, m, m, size1, i, 1, m);
            block_mult_block(h, y, z, size, m, size1, m);
            get_block(aa, x, m, size, size1, s, 1, m);
            block_minus_block(x, z, y, size, size1, m);
            put_block(aa, y, m, size, size1, s, 1, m);
        }
    }

    for (j = 1; j <= max_cols; j++)
    {
        global_j = (j - 1) * p + q + 1;

        if (j == max_cols && l > 0 && global_j == k) size1 = l;
        else size1 = m;

        for (s = i + 1; s <= k; ++s)
        {
            if (s != k) size = m;
            else if (l > 0) size = l;
            else size = m;

            get_block(buf, h, m, size, m, s, 1, m);
            get_block(b, y, m * max_cols, m, size1, i, j, m);
            block_mult_block(h, y, z, size, m, size1, m);
            get_block(b, x, m * max_cols, size, size1, s, j, m);
            block_minus_block(x, z, y, size, size1, m);
            put_block(b, y, m * max_cols, size, size1, s, j, m);
        }
    }
}

/* обратный ход простого варианта */
inline void answer1(double *a, double *e, int n, int m)
{
    int i, j, k;

    for (i = n - 2; i >= 0; --i)
        for (j = n - 1; j > i; --j)
        {
            if ((i != 0) && (j == n - 1))
            {
                for (k = 0; k <= i; ++k)
                        e[i * m + k] -= a[i * m + j] * e[j * m + k];
                for (k = i + 1; k < n; ++k)
                    e[i * m + k] = -a[i * m + j] * e[j * m + k];
            }
            else
            {
                for (k = 0; k < n; ++k)
                    e[i * m + k] -= a[i * m + j] * e[j * m + k];
            }
            a[i * m + j] = 0;
        }
}

/* обратный ход MPI варианта */
void answer(double *a, double *b, double *x, double *y, double *z, double *h, double *buf, int n, int m, int l, int k, int max_cols, int j)
{
    int q, p, size, size1;
    int i, t, local_j, global_t;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (j == k && l > 0) size = l;
    else size = m;

    if ((j - 1) % p == q)
    {
        local_j = (j - 1) / p + 1;

        get_block(a + (local_j - 1) * m * n, buf, m, n, m, 1, 1, m);

        MPI_Bcast(buf, n * m, MPI_DOUBLE, q, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast(buf, n * m, MPI_DOUBLE, (j - 1) % p, MPI_COMM_WORLD);
    }

    for (i = j - 1; i > 0; --i)
    {
        get_block(buf, h, m, m, size, i, 1, m);

        for (t = 1; t <= max_cols; t++)
        {
            global_t = (t - 1) * p + q + 1;

            if (t == max_cols && l > 0 && global_t == k) size1 = l;
            else size1 = m;

            get_block(b, y, m * max_cols, size, size1, j, t, m);
            block_mult_block(h, y, z, m, size, size1, m);
            get_block(b, x, m * max_cols, m, size1, i, t, m);
            block_minus_block(x, z, y, m, size1, m);
            put_block(b, y, m * max_cols, m, size1, i, t, m);
        }
    }
}

/* норма матрицы */
double norm_matrix(double *a, int n, int m, int max_cols)
{
    int q, p;
    double res, max;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    res = norm_matrix_q(a, n, m, max_cols);
    MPI_Allreduce(&res, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return max;
}

/* вклад в норму матрицы от процесса с номером q */
double norm_matrix_q(double *a, int n, int m, int max_cols)
{
    int i, j;
    double s = 0, max = 0;

    for (j = 0; j < m; j++)
    {
        for (i = 0; i < n * max_cols; ++i)
        {
            s += fabs(a[i * m + j]);
            if ((i + 1) % n == 0 && i > 0)
            {
                if (s > max) max = s;
                s = 0;
            }
        }
    }
    return max;
}

int solve(double *a, double *b, double *x, double *y, double *z, double *h, double *buf, int n, int m, int k, int l, int max_cols, double z1)
{
    int q, p, size, size1;
    int i, j, local_i, global_j, ff, local_ff = 1;
    double *aa;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* где цикл до k - там нужен пересчет в локальных индексах;
       где до max_cols - не нужен пересчет; */

    /* прямой ход */
    for (i = 1; i <= k; ++i)
    {
        if (i != k) size = m;
        else if (l > 0) size = l;
        else size = m;

        if ((i - 1) % p == q) /* отправить (i,i)-тый блок остальным */
        {
            local_i = (i - 1) / p + 1;

            get_block(a + (local_i - 1) * m * n, buf, m, n, m, 1, 1, m);

            MPI_Bcast(buf, n * m, MPI_DOUBLE, q, MPI_COMM_WORLD);
        }
        else /* принять блок */
        {
            MPI_Bcast(buf, n * m, MPI_DOUBLE, (i - 1) % p, MPI_COMM_WORLD);
        }

        get_block(buf, x, m, size, size, i, 1, m);

        put_block_E(y, m, m);

        /* найти обратную к x и положить в y */
        int res = solve1(x, y, size, m, z1);

        if (res < 0)
        {
            return -1;
        }
        else
        {
            ff = f(i, q + 1, p);
            local_ff = (ff - 1) / p + 1;

            /* умножаем на обратную i-тую строку матрицы A*/
            for (j = local_ff; j <= max_cols; j++)
            {
                global_j = (j - 1) * p + q + 1;

                if (j == max_cols && l > 0 && global_j == k) size1 = l;
                else size1 = m;

                aa = a + (j - 1) * n * m;

                get_block(aa, x, m, size, size1, i, 1, m);
                block_mult_block(y, x, z, size, size, size1, m);
                put_block(aa, z, m, size, size1, i, 1, m);
            }

            /* умножаем на обратную i-тую строку присоединенной матрицы*/
            for (j = 1; j <= max_cols; j++)
            {
                global_j = (j - 1) * p + q + 1;

                if (j == max_cols && l > 0 && global_j == k) size1 = l;
                else size1 = m;

                get_block(b, x, m * max_cols, size, size1, i, j, m);
                block_mult_block(y, x, z, size, size, size1, m);
                put_block(b, z, m * max_cols, size, size1, i, j, m);
            }
        }

        /* цикл вычитания */
        if (i < k) subtraction(a, b, x, y, z, h, buf, n, m, l, k, max_cols, i, local_ff);
    }

    /* обратный ход */
    for (j = k; j > 1; --j)
    {
        answer(a, b, x, y, z, h, buf, n, m, l, k, max_cols, j);
    }
    return 0;
}

int solve1(double *a, double *e, int n, int m, double z1)
{
    int i, res;
    double z = norm_matrix1(a,n,m);

    for (i = 0; i < n; ++i)
    {
      res = division1(a, e, i, n, z, m, z1);
      if (res < 0) return -1;
      if (i != n - 1) subtraction1(a, e, i, n, m);
    }
    answer1(a, e, n, m);
    return 0;
}

inline double norm_matrix1(double *a, int n, int m)
{
    int i, j;
    double s = 0, max = 0;

    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
            s += fabs(a[i * m + j]);
        if (s > max) max = s;
        s = 0;
    }
    return max;
}

void function(void *a, void *b, int *len, MPI_Datatype *datatype)
{
    int i;
    (void) datatype;
    double *aa = (double*) a;
    double *bb = (double*) b;
    for (i = 0; i < *len; i++) bb[i] = aa[i] + bb[i];
}

double residual(double *a, double *b, double *r, double *x, double *y, double *z, double *h, int n, int m, int k, int l, int max_cols)
{
    int q, p, size, size1, size2;
    int i, j, t, v, w, g;
    int local_j, local_t;
    double s, max = 0;
    MPI_Op MPI_MATRIX_PLUS_MATRIX;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    matrix H(n,m);

    for (j = 1; j <= k; j++)
    {
        if (j == k && l > 0) size = l;
        else size = m;

        if ((j - 1) % p == q)
        {
            local_j = (j - 1) / p + 1;

            get_block(b, H.get_array(), m * max_cols, n, m, 1, local_j, m);

            /* отослать по блокам, свои умножить */
            MPI_Bcast (H.get_array(), n * m, MPI_DOUBLE, q, MPI_COMM_WORLD);
        }
        else
        {
            /* принять по блокам, свои умножить */
            MPI_Bcast (H.get_array(), n * m, MPI_DOUBLE, (j - 1) % p, MPI_COMM_WORLD);
        }

        // теперь у всех в H.get_array() лежит j-тый столбец матрицы B

        for (i = 1; i <= k; i++)
        {
            if (i == k && l > 0) size1 = l;
            else size1 = m;

            for (t = q; t < k; t += p)
            {
                if (t == k - 1 && l > 0) size2 = l;
                else size2 = m;

                local_t = t / p + 1;

                get_block(a + (local_t - 1) * n * m, x, m, size1, size2, i, 1, m);
                get_block(H.get_array(), y, m, size2, size, t + 1, 1, m);
                block_mult_block(x, y, z, size1, size2, size, m);
                get_block(r, h, m, size1, m, i, 1, m);
                block_plus_block(h, z, x, size1, m, m);
                put_block(r, x, m, size1, m, i, 1, m);
            }
        }

        // теперь у каждого в r блоки, сложив которые получим j-тый блочный столбец матрицы AxB

        MPI_Op_create((MPI_User_function*) function, 1, &MPI_MATRIX_PLUS_MATRIX);
        MPI_Allreduce(r, H.get_array(), n * m, MPI_DOUBLE, MPI_MATRIX_PLUS_MATRIX, MPI_COMM_WORLD);
        MPI_Op_free(&MPI_MATRIX_PLUS_MATRIX);

        for (v = 0; v < size; v++)
        {
            s = 0;
            for (w = 0; w < n; w++)
            {
                s += H.get_array()[w * m + v];
            }
            s--; // вычесть единичную матрицу
            if (j == 1 && v == 0) max = s;
            else if (s > max) max = s;
        }

        for (g = 0; g < n * m; g++) r[g] = 0;
    }
    return max;
}
