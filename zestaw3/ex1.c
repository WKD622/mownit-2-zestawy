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

void print_matrix(int size, double **matrix, char *title)
{
    int i, j;
    printf("%s\n", title);
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

double **matrix_1d_to_2d(int size, gsl_matrix *matrix)
{
    int i, j, k = 0;
    double **matrix_2d = allocate_matrix(size);
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            matrix_2d[i][j] = matrix->data[k];
            k++;
        }
    }
    return matrix_2d;
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

void fill_vector_with_zeros(int size, double *vector)
{
    int i;
    for (i = 0; i < size; i++)
        vector[i] = 0;
}

float arithmetic_average(int size, double *x)
{
    double sum = 0.0;
    int i;
    for (i = 0; i < size; i++)
        sum = sum + x[i];
    if (sum != 0.0)
        return sum / size;
    else
        return 0.0;
}

void jacobian_method(int size, double **A, gsl_vector *b, double accuracy)
{
    int i, j, s, k;
    double **L_U, **D;

    // TWORZE MACIERZE A = L + D + U
    L_U = allocate_matrix(size);
    D = allocate_matrix(size);
    fill_matrix_with_zeros(size, L_U);
    fill_matrix_with_zeros(size, D);

    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            if (i == j)
            {
                D[i][j] = A[i][j];
            }
            else
            {
                L_U[i][j] = A[i][j];
            }
        }
    }
    print_matrix(size, L_U, "L+U");
    // ODWRACAM MACIERZ D
    double *D_1d = matrix_2d_to_1d(size, D);
    double **D_after_inversion;
    double inva[size * size];

    gsl_matrix_view m = gsl_matrix_view_array(D_1d, size, size);
    gsl_matrix_view inv = gsl_matrix_view_array(inva, size, size);
    gsl_permutation *p = gsl_permutation_alloc(size);

    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_LU_invert(&m.matrix, p, &inv.matrix);

    D_after_inversion = matrix_1d_to_2d(size, &inv.matrix);

    print_matrix(size, D_after_inversion, "D after inversion");

    for (i = 0; i < size * size; i++)
    {
        if (inv.matrix.data[i] != 0.0)
            inv.matrix.data[i] = inv.matrix.data[i] * (-1);
    }

    D_after_inversion = matrix_1d_to_2d(size, &inv.matrix);
    print_matrix(size, D_after_inversion, "D after inversion * -1");

    // WYLICZAM MACIERZ M = L_U * D^(-1)
    double *L_U_1d = matrix_2d_to_1d(size, L_U);
    double **M_2d = allocate_matrix(size);
    gsl_matrix_view L_U_v = gsl_matrix_view_array(L_U_1d, size, size);
    gsl_matrix *M = gsl_matrix_alloc((size_t)size, (size_t)size);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &inv.matrix, &L_U_v.matrix, 0.0, M);
    M_2d = matrix_1d_to_2d(size, M);
    print_matrix(size, M_2d, "D^(-1) * L_U");

    // TERAZ MAM JUZ MACIERZ M ORAZ b, PRZYBLIZAM ROZWIAZANIE
    double *x = allocate_vector(size);
    double *x_copy = allocate_vector(size);
    double *x_dif = allocate_vector(size);
    fill_vector_with_zeros(size, x);
    fill_vector_with_zeros(size, x_dif);
    int number_of_iterations = 0;

    while (1)
    {
        for (j = 0; j < size; j++)
        {
            x_copy[j] = x[j];
            x_dif[j] = x[j];
        }
        // OD X1 DO XN
        for (j = 0; j < size; j++)
        {
            x[j] = b->data[j] / A[j][j];
            for (k = 0; k < size; k++)
            {
                if (k != j)
                    x[j] = x[j] + M_2d[j][k] * x_copy[k];
            }
        }
        number_of_iterations++;
        for (j = 0; j < size; j++)
        {
            x_dif[j] = x_dif[j] - x[j];
            if (x_dif[j] < 0.0)
                x_dif[j] = x_dif[j] * (-1.0);
        }
        if (arithmetic_average(size, x_dif) < accuracy)
        {
            printf("Number of iterations %d\n\n", number_of_iterations);
            for (i = 0; i < size; i++)
                printf("x%d = %f\n", i, x[i]);

            free_matrix(size, L_U);
            free_matrix(size, D);
            free_matrix(size, D_after_inversion);
            free_matrix(size, M_2d);
            free_vetor(D_1d);
            free_vetor(L_U_1d);
            free_vetor(x);
            free_vetor(x_dif);
            free(p);
            free(M);
            return;
        }
    }
}

int main()
{
    srand(time(NULL));
    int size = 10;
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

    // gsl_vector_fprintf(stdout, b, "%g");

    print_matrix(size, A, "A");

    // JACOBI
    jacobian_method(size, A, b, 0.0000000001);

    free_vetor(x);
    free_vetor(A_1d);
    free_matrix(size, A);
    free(b);

    return 0;
}