#pragma once
#include <memory>
#include <stdio.h>
#include <stdlib.h>

class matrix
{
private:
    double *a;
    int size;
public:
    matrix() { a = 0; size = 0; }
    matrix(int n)
    {
        size = n;
        a = new double[size * size];

        for (int i = 0; i < size * size; i++) a[i] = 0;
    }
    matrix(int n, int p)
    {
        size = n;
        a = new double[size * p];

        for (int i = 0; i < size * p; i++) a[i] = 0;
    }
    matrix(int n, int m, int p)
    {
        size = m * p;
        a = new double[n * m * p];

        for (int i = 0; i < n * m * p; i++) a[i] = 0;
    }
    ~matrix() { if (a) delete [] a; size = 0;}
    double* get_array() { return a; }
    double get_ij(int i, int j) {return a[i * size + j];}
    double get(int i)const {return a[i];}
    void set_ij(int i, int j, double x) {a[i * size + j] = x;}
    void set(int i, double x) {a[i] = x;}
    int get_size()const { return size; }
    int read(char *filename);
    void print();
    void print(int m, int l);
    void init();
    void init_diag();
    double f(int i, int j, int n);
};

class array
{
private:
    double *a;
    int size;
public:
    array() { a = 0; size = 0; }
    array(int n)
    {
        size = n;
        a = new double[size];
    }
    ~array() { delete [] a; size = 0;}
    double* get_array() { return a; }
    double get(int i)const {return a[i];}
    void set(int i, double x) {a[i] = x;}
    int get_size()const { return size; }
    void print()
    {
        for (int i = 0; i < size; i++) printf("%.3f ", a[i]);
        printf("\n");
    }
    void init()
    {
        for (int i = 0; i < size; i++) a[i] = 0;
    }
};

