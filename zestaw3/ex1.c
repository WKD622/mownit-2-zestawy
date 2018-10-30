#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_linalg.h>

double **allocate_matrix(int size)
{
    int i;
    double **matrix = (double **)malloc(size * sizeof(double *));
    for (i = 0; i < size; i++)
    {
        matrix[i] = (double *)malloc(size * sizeof(double));
    }
    return matrix;
}

void fill_matrix(int size, double **matrix)
{

    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            if (i == j)
            {
                if (i == 0 || i == size - 1)
                    matrix[i][j] = 1;
                else
                    matrix[i][j] = 2;
            }
            else if (j == i + 1)
            {
                matrix[i][j] = 1 / ((double)j + 1.0);
            }
            else if (i == j + 1)
            {
                matrix[i][j] = 1 / ((double)i + 1.0);
            }
            else
                matrix[i][j] = 0;
        }
    }
}

void print_matrix(int size, double **matrix)
{
    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            printf("%12f", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void free_matrix(int size, double **matrix)
{
    int i;
    for (i = 0; i < size; i++)
        free(matrix[i]);
    free(matrix);
}

double *allocate_vector(int size)
{
    double *vector = (double *)malloc(size * sizeof(double));
    return vector;
}

void fill_vector(int size, double *vector)
{
    int i;
    for (i = 0; i < size; i++)
    {
        vector[i] = rand() % 2;
    }
}

void print_vector(int size, double *vector)
{
    int i;
    for (i = 0; i < size; i++)
        printf("%12f", vector[i]);
    printf("\n");
}

void free_vetor(double *vector)
{
    free(vector);
}

double *matrix_2d_to_1d(int size, double **matrix)
{
    double *matrix_1d = allocate_vector(size * size);
    int i, j, k = 0;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            matrix_1d[k] = matrix[i][j];
            k++;
        }
    }
    return matrix_1d;
}

void fill_matrix_with_zeros(int size, double **matrix)
{
    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            matrix[i][j] = 0;
        }
    }
}

void jacobian_method(int size, gsl_matrix_view A, gsl_vector_view b)
{
    int i, j;
    // gsl_matrix L, D, U;
    // for (i = 0; i < size; i++)
    // {
    //     for (j = 0; j < size; j++)
    //     {
    //         if (i = j)
    //         {
    //         }
    //     }
    // }
}

int main()
{
    srand(time(NULL));
    int size = 5;
    double **A = allocate_matrix(size);
    double *x = allocate_vector(size);
    gsl_vector *b = gsl_vector_alloc((size_t)size);

    fill_matrix(size, A);
    fill_vector(size, x);
    // print_matrix(size, A);
    // print_vector(size, x);

    double *A_1d = matrix_2d_to_1d(size, A);
    // print_vector(size * size, A_1d);

    gsl_matrix_view A_v = gsl_matrix_view_array(A_1d, (size_t)size, (size_t)size);
    gsl_vector_view x_v = gsl_vector_view_array(x, (size_t)size);

    // wyliczanie wektora b
    gsl_blas_dgemv(CblasNoTrans, 1.0, &A_v.matrix, &x_v.vector, 0.0, b);

    gsl_vector_fprintf(stdout, b, "%g");

    int i, j;
    for (i = 0; i < size * size; i++)
    {
        printf("%12f", A_v.matrix.data[i]);
    }
    free_vetor(x);
    free_vetor(A_1d);
    free_matrix(size, A);
    free(b);

    return 0;
}