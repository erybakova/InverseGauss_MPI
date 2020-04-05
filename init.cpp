#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/sysinfo.h>

#include "mpi.h"
#include "matrix.h"
#include "solve.h"

#define MAX_PRINT 10

void print_matrix_0(double *a, int n, int m, int l, int k)
{
    int p, size, t, N, K;
    int i, j, h, max_blocks_num, max_num_block_size;
    double *buf;

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Status status;

    for (max_blocks_num = 0; ; max_blocks_num++) if (max_blocks_num * m > MAX_PRINT) break;
    max_num_block_size = MAX_PRINT - (max_blocks_num - 1) * m;

    N = (n > MAX_PRINT ? MAX_PRINT : n);
    K = (k < max_blocks_num ? k : max_blocks_num);

    //printf("will stop printing at %d block with size %d, K = %d\n", max_blocks_num, max_num_block_size, K);

    buf = new double[m];

    for (i = 0; i < N; i++)
    {
        h = 0;

        for (j = 0; j < K; j++)
        {
            if (j == k - 1)
            {
                if (j == max_blocks_num - 1)
                {
                    if (l > 0) size = (l < max_num_block_size ? l : max_num_block_size);
                    else size = max_num_block_size;
                }
                else if (l > 0) size = l;
                else size = m;
            }
            else if (j == max_blocks_num - 1) size = max_num_block_size;
            else size = m;

            if (j % p > 0)
            {
                MPI_Recv(buf, m, MPI_DOUBLE, j % p, 0, MPI_COMM_WORLD, &status);

                for (t = 0; t < size; t++) printf("%12.6f ", buf[t]);
            }
            else
            {
                for (t = 0; t < size; t++)
                {
                    printf("%12.6f ", (a + h * m + i * m)[t]);
                }

                h += n;
            }
        }
        printf("\n");
    }
    printf("\n");

    delete [] buf;
}

void print_matrix_others(double *a, int n, int m, int l, int k)
{
    int q, p, size, max_blocks_num, max_num_block_size, N, K;
    int i, j, h;
    double *buf;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    for (max_blocks_num = 0; ; max_blocks_num++) if (max_blocks_num * m > MAX_PRINT) break;

    max_num_block_size = MAX_PRINT - (max_blocks_num - 1) * m;

    N = (n > MAX_PRINT ? MAX_PRINT : n);
    K = (k < max_blocks_num ? k : max_blocks_num);

    buf = new double[m];

    for (i = 0; i < N; i++)
    {
        h = 0;

        for (j = q; j < K; j += p)
        {
            if (j == k - 1)
            {
                if (j == max_blocks_num - 1)
                {
                    if (l > 0) size = (l < max_num_block_size ? l : max_num_block_size);
                    else size = max_num_block_size;
                }
                else if (l > 0) size = l;
                else size = m;
            }
            else if (j == max_blocks_num - 1) size = max_num_block_size;
            else size = m;

            MPI_Send(a + h * m + i * m, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            h += n;
        }
    }

    delete [] buf;
}

void print_matrix(double *a, int n, int m, int l, int k)
{
    int q;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);

    if (q > 0) print_matrix_others(a, n, m, l, k);
    else print_matrix_0(a, n, m, l, k);
}

void print_one(double *a, int n, int m, int max_cols)
{
    int i, j;

    for (i = 0; i < n * max_cols; i++)
    {
        for (j = 0; j < m; j++)
            printf("%.5f ", a[i * m + j]);
        printf("\n");
    }
    printf("\n");
}

int read_string(FILE *fp, double *buf, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (fscanf(fp, "%lf", &buf[i]) != 1) return 0;
    }

    return 1;
}

int init_from_file_0(double *a, int n, int m, int l, int k, char *filename)
{
    int p, res, size;
    int i, j, h;
    double *buf;
    FILE *fp;

    MPI_Comm_size(MPI_COMM_WORLD, &p);

    fp = fopen(filename, "r");

    MPI_Bcast(&fp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!fp)
    {
        return -1;
    }

    buf = new double[n];

    for (i = 0; i < n; i++)
    {
        h = 0;

        res = read_string(fp, buf, n);

        MPI_Bcast(&res, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (res == 0)
        {
            delete [] buf;
            return -2;
        }

        for (j = 0; j < k; j++)
        {
            if ((j == k - 1) && (l > 0)) size = l;
            else size = m;

            if (j % p > 0)
            {
                MPI_Send(buf + j * m, size, MPI_DOUBLE, j % p, 0, MPI_COMM_WORLD);
            }
            else
            {
                memcpy(a + m * h + m * i, buf + j * m, size * sizeof(double));
                h += n;
            }
         }
    }

    delete [] buf;

    return 1;
}

int init_from_file_others(double *a, int n, int m, int l, int k)
{
    int q, p, size;
    int i, j, fp, h;
    int res;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Bcast(&fp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!fp) return -1;

    for (i = 0; i < n; i++)
    {
        h = 0;

        MPI_Bcast(&res, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!res) return -2;

        for (j = q; j < k; j += p)
        {
            if ((j == k - 1) && (l > 0)) size = l;
            else size = m;

            MPI_Recv(a + m * h + m * i, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

            h += n;
        }
    }

    return 1;
}

int init_from_file(double *a, int n, int m, int l, int k, char *filename)
{
    int q;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);

    if (q == 0) return init_from_file_0(a, n, m, l, k, filename);
    else return init_from_file_others(a, n, m, l, k);
}

double f_(int i, int j, int n)
{
    //return fabs(i - j);
    //return 1. /(i + j + 1);
    int k = i > j ? i : j;
    return n - k;
}

void init_by_formula(double *a, int n, int m, int l, int k)
{
    int q, p, size;
    int i, j, t, h;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    for (i = 0; i < n; i++)
    {
        h = 0;

        for (j = q; j < k; j += p)
        {
            if ((j == k - 1) && (l > 0)) size = l;
            else size = m;

            for (t = 0; t < size; t++)
            {
                (a + h * m + i * m)[t] = f_(i, j * m + t, n);
            }

            h += n;
        }
    }
}

void init_diag(double *a, int n, int m, int l, int k)
{
    int q, p, size, max_cols;
    int i, j, t, h;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    max_cols = (k % p > 0 ? k / p + 1 : k / p);

    for (i = 0; i < n; i++)
    {
        h = 0;

        for (j = q; j < k; j += p)
        {
            if ((j == k - 1) && (l > 0)) size = l;
            else size = m;

            for (t = 0; t < size; t++)
            {
                if (j * m + t == i) (a + i * m * max_cols + h)[t] = 1;
            }

            h += m;
        }
    }
}

void print_matrix_0_1(double *a, int n, int m, int l, int k)
{
    int p, size, t, max_cols, N, K;
    int i, j, h, max_blocks_num, max_num_block_size;
    double *buf;

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Status status;

    for (max_blocks_num = 0; ; max_blocks_num++) if (max_blocks_num * m > MAX_PRINT) break;
    max_num_block_size = MAX_PRINT - (max_blocks_num - 1) * m;
    max_cols = (k % p > 0 ? k / p + 1 : k / p);

    N = (n > MAX_PRINT ? MAX_PRINT : n);
    K = (k < max_blocks_num ? k : max_blocks_num);

    //printf("will stop printing at %d block with size %d, K = %d\n", max_blocks_num, max_num_block_size, K);

    buf = new double[m];

    for (i = 0; i < N; i++)
    {
        h = 0;

        for (j = 0; j < K; j++)
        {
            if (j == k - 1)
            {
                if (j == max_blocks_num - 1)
                {
                    if (l > 0) size = (l < max_num_block_size ? l : max_num_block_size);
                    else size = max_num_block_size;
                }
                else if (l > 0) size = l;
                else size = m;
            }
            else if (j == max_blocks_num - 1) size = max_num_block_size;
            else size = m;

            if (j % p > 0)
            {
                MPI_Recv(buf, m, MPI_DOUBLE, j % p, 0, MPI_COMM_WORLD, &status);

                for (t = 0; t < size; t++) printf("%12.6f ", buf[t]);
            }
            else
            {
                for (t = 0; t < size; t++)
                {
                    printf("%12.6f ", (a + i * m * max_cols + h)[t]);
                }

                h += m;
            }
        }
        printf("\n");
    }
    printf("\n");

    delete [] buf;
}

void print_matrix_others_1(double *a, int n, int m, int l, int k)
{
    int q, p, size, max_cols, max_blocks_num, max_num_block_size, N, K;
    int i, j, h;
    double *buf;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    for (max_blocks_num = 0; ; max_blocks_num++) if (max_blocks_num * m > MAX_PRINT) break;

    max_num_block_size = MAX_PRINT - (max_blocks_num - 1) * m;

    N = (n > MAX_PRINT ? MAX_PRINT : n);
    K = (k < max_blocks_num ? k : max_blocks_num);

    max_cols = (k % p > 0 ? k / p + 1 : k / p);

    buf = new double[m];

    for (i = 0; i < N; i++)
    {
        h = 0;

        for (j = q; j < K; j += p)
        {
            if (j == k - 1)
            {
                if (j == max_blocks_num - 1)
                {
                    if (l > 0) size = (l < max_num_block_size ? l : max_num_block_size);
                    else size = max_num_block_size;
                }
                else if (l > 0) size = l;
                else size = m;
            }
            else if (j == max_blocks_num - 1) size = max_num_block_size;
            else size = m;

            MPI_Send(a + i * m * max_cols + h, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

            h += m;
        }
    }

    delete [] buf;
}

void print_matrix_1(double *a, int n, int m, int l, int k)
{
    int q;

    MPI_Comm_rank(MPI_COMM_WORLD, &q);

    if (q > 0) print_matrix_others_1(a, n, m, l, k);
    else print_matrix_0_1(a, n, m, l, k);
}
