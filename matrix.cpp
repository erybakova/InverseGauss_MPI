#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include "math.h"
#include "matrix.h"

#define MAX_PRINT 10

int matrix:: read(char *filename)
{
    FILE *fp;
    if (!(fp = fopen(filename, "r"))) return -1;
    for (int i = 0; i < size * size; i++)
        {if (fscanf(fp, "%lf", a + i) != 1)
          {fclose(fp); return -2;}}
    fclose(fp);
    return 0;
}

void matrix:: print()
{
    int mp = (size > MAX_PRINT ? MAX_PRINT : size);
    for (int i = 0; i < mp; i++)
    {
        for (int j = 0; j < mp; j++)
           {printf(" %.3f", *(a + i * size + j));}
        printf("\n");
    }
    printf("\n");
}

void matrix:: print(int m, int l)
{
    int mp = (m > MAX_PRINT ? MAX_PRINT : m);
    int lp = (l > MAX_PRINT ? MAX_PRINT : l);
    for (int i = 0; i < mp; i++)
    {
        for (int j = 0; j < lp; j++)
           {printf(" %.3f", *(a + i * l + j));}
        printf("\n");
    }
    printf("\n");
}

void matrix:: init()
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
              {a[i * size + j] = f(i, j, size);}
}

void matrix:: init_diag()
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
        {
            if (j == i) a[i * size + j] = 1;
            else a[i * size + j] = 0;
        }
}

double matrix:: f(int i, int j, int n)
{
    //return fabs(i - j);
    //return 1. /(i + j + 1);
    int k = i > j ? i : j;
    return n - k;
}
