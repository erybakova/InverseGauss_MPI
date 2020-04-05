#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
#include "init.h"
#include "matrix.h"
#include "solve.h"

int main(int argc, char **argv)
{
    int n, m, q, p;
    int k, l, max_cols;
    int res;
    char *filename;
    double t, r, z1;

    if ((argc < 3) || (argc > 4) || ((n = atoi(argv[1])) <= 0)
        || ((m = atoi(argv[2])) <= 0) || ((m = atoi(argv[2])) > n))
    {
        printf("Usage: %s n m <filename>\n", argv[0]);
        return 1;
    }

    filename = argv[3];

    MPI_Init(&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &q); /* номер текущего процесса в группе всех процессов */
    MPI_Comm_size (MPI_COMM_WORLD, &p); /* общее количество запущенных процессов */

    k = (n % m > 0 ? n / m + 1 : n / m); /* число блочных строк матрицы */
    max_cols = (k % p > 0 ? k / p + 1 : k / p); /* максимальное количество блочных столбцов на процесс */
    l = n - (n / m) * m; /* остаточек */

    //if (q == 0) printf("n, m, l: %d, %d, %d;\nk and p: %d, %d;\nmax_cols = %d\n", n, m, l, k, p, max_cols);

    matrix A(n, m, max_cols);
    matrix B(n, m, max_cols);
    matrix X(m);
    matrix Y(m);
    matrix Z(m);
    matrix H(m);
    matrix T(n, m);
    matrix R(n, m);

    if (filename)
    {
        res = init_from_file(A.get_array(), n, m, l, k, filename);
        if (res < 0)
        {
            switch (res)
            {
                case -1: {if (q == 0) printf("Cannot open %s!\n", filename); break;}
                case -2: {if (q == 0) printf("Cannot read file %s!\n", filename); break;}
                default: {if (q == 0) printf("Unknown error %d in %s!\n", res, filename); break;}
            }
            MPI_Finalize();
            return 0;
        }
    }
    else init_by_formula(A.get_array(), n, m, l, k);

    init_diag(B.get_array(), n, m, l, k);

    /*printf("\nPart of matrix A for %ds process:\n", q);
    print_one(A.get_array(), n, m, max_cols);*/

    if (q == 0) printf("\nMatrix A:\n");
    print_matrix(A.get_array(), n, m, l, k);

    z1 = norm_matrix(A.get_array(), n, m, max_cols);

    t = MPI_Wtime();

    res = solve(A.get_array(), B.get_array(), X.get_array(), Y.get_array(), Z.get_array(), H.get_array(), T.get_array(), n, m, k, l, max_cols, z1);

    MPI_Barrier(MPI_COMM_WORLD);

    t = (MPI_Wtime() - t);

    if (res < 0)
    {
        if (q == 0) printf("Method is unworkable!\n");
        MPI_Finalize();
        return 0;
    }
    else
    {
        if (q == 0) printf("\nMatrix A^(-1):\n");
        print_matrix_1(B.get_array(), n, m, l, k);
    }

    if ((n < 10000) || (p > 1))
    {
        if (filename)
        {
            res = init_from_file(A.get_array(), n, m, l, k, filename);
        }
        else init_by_formula(A.get_array(), n, m, l, k);

        r = residual(A.get_array(), B.get_array(), R.get_array(), X.get_array(), Y.get_array(), Z.get_array(), H.get_array(), n, m, k, l, max_cols);

        if (q == 0) printf("Residual = %.e, n = %d, m = %d, p = %d, elapsed: %.2f\n", r, n, m, p, t);
    }
    else if (q == 0) printf("Residual = %.e, n = %d, m = %d, p = %d, elapsed: %.2f\n", 0e0, n, m, p, t);

    MPI_Finalize();

    return 0;
}
