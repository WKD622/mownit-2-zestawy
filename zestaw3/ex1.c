#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <string.h>

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

void identity_matrix(int size, double **M)
{
    int i, j;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            if (i == j)
                M[i][j] = 1;
            else
                M[i][j] = 0;
        }
    }
}

gsl_vector *copy_one_vector_to_another(int size, gsl_vector *a)
{
    gsl_vector *b = gsl_vector_alloc(size);
    int i;
    for (i = 0; i < size; i++)
        b->data[i] = a->data[i];
    return b;
}

void jacobian_method(int size, double **A, gsl_vector *b, double accuracy)
{
    int i, j, s, k;
    double **L_U, **D;

    // print_matrix(size, A, "A");

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
    // print_matrix(size, L_U, "L+U");
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

    // print_matrix(size, D_after_inversion, "D after inversion");

    for (i = 0; i < size * size; i++)
    {
        if (inv.matrix.data[i] != 0.0)
            inv.matrix.data[i] = inv.matrix.data[i] * (-1);
    }

    D_after_inversion = matrix_1d_to_2d(size, &inv.matrix);
    // print_matrix(size, D_after_inversion, "D after inversion * -1");

    // WYLICZAM MACIERZ M = L_U * D^(-1)
    double *L_U_1d = matrix_2d_to_1d(size, L_U);
    double **M_2d = allocate_matrix(size);
    gsl_matrix_view L_U_v = gsl_matrix_view_array(L_U_1d, size, size);
    gsl_matrix *M = gsl_matrix_alloc((size_t)size, (size_t)size);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &inv.matrix, &L_U_v.matrix, 0.0, M);
    M_2d = matrix_1d_to_2d(size, M);
    // print_matrix(size, M_2d, "D^(-1) * L_U");

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
            free_vetor(x_copy);
            free(p);
            free(M);
            return;
        }
    }
}

void chebyshev_method(int size, double **A, gsl_vector *b, double *x0, int iter_num, double l_max, float l_min)
{
    int i, j, s;
    float alpha, beta;
    float d = (l_max + l_min) / 2.0;
    float c = (l_max - l_min) / 2.0;
    double **pre_cond = allocate_matrix(size);
    double *x_1 = allocate_vector(size);
    double *A_1d = matrix_2d_to_1d(size, A);
    gsl_permutation *pom = gsl_permutation_alloc(size);
    identity_matrix(size, pre_cond);
    double *pre_cond_1d = matrix_2d_to_1d(size, pre_cond);
    gsl_matrix_view pre_cond_1d_m = gsl_matrix_view_array(pre_cond_1d, size, size);
    gsl_vector *r = copy_one_vector_to_another(size, b);
    gsl_vector *z = gsl_vector_alloc(size);
    gsl_vector *p;
    gsl_vector *result = gsl_vector_alloc((size_t)size);
    double *A_copy = allocate_vector(size * size);

    for (i = 0; i < size; i++)
        x_1[i] = x0[i];

    gsl_vector_view x = gsl_vector_view_array(x_1, size);

    for (i = 0; i < 1000; i++)
    {
        memcpy(A_copy, A_1d, size * size * sizeof(double));
        gsl_matrix_view A_view = gsl_matrix_view_array(A_copy, size, size);

        gsl_linalg_LU_decomp(&pre_cond_1d_m.matrix, pom, &s);
        gsl_linalg_LU_solve(&pre_cond_1d_m.matrix, pom, r, z);

        if (i == 0)
        {
            p = copy_one_vector_to_another(size, z);
            alpha = 1.0 / d;
        }
        else
        {
            beta = pow(c * alpha / 2.0, 2.0);
            alpha = 1.0 / (d - beta / alpha);
            for (j = 0; j < size; j++)
                p->data[j] = p->data[j] * beta + z->data[j];
        }
        for (j = 0; j < size; j++)
            x.vector.data[j] = x.vector.data[j] + beta * p->data[j];

        gsl_blas_dgemv(CblasNoTrans, 1.0, &A_view.matrix, &x.vector, 0.0, result);
        for (j = 0; j < size; j++)
            r->data[j] = b->data[j] - result->data[j];
    }
    print_vector(size, x_1);

    free_matrix(size, pre_cond);
    free_vetor(x_1);
    free(pom);
    free(z);
    free(r);
    free_vetor(pre_cond_1d);
    free_vetor(A_copy);
    free_vetor(A_1d);

    // function [x] =  SolChebyshev002(A,b,x0,iterNum,lMax,lMin)

    //   d=(lMax+lMin)/2;
    //   c=(lMax-lMin)/2;
    //   preCond=eye(size(A)); %preconditioner
    //   x=x0;
    //   r=b-A*x;

    //   for i = 1:iterNum % size(A,1)
    //       z = linsolve(preCond,r);
    //       if (i==1)
    //           p=z;
    //           alpha=1/d;
    //       else
    //           beta=(c*alpha/2)^2;
    //           alpha=1/(d - beta/alpha);
    //           p=z+beta*p;
    //       end;

    //       x=x+alpha*p;
    //       r=b-A*x; %(=r-alpha*A*p)
    //       if (norm(r)<1e-15), break; end; %stop if necessary [Given A on the order of eps this will make the iteration stop too soon. Should be norm(r)<eps*norm(b) or norm(r)<eps*norm(A*x) depending on the details that I am unaware of; i.e., "1e-15" carries units.]
    //   end;
    // end
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

    double *A_1d = matrix_2d_to_1d(size, A);

    gsl_matrix_view A_v = gsl_matrix_view_array(A_1d, (size_t)size, (size_t)size);
    gsl_vector_view x_v = gsl_vector_view_array(x, (size_t)size);

    // wyliczanie wektora b
    gsl_blas_dgemv(CblasNoTrans, 1.0, &A_v.matrix, &x_v.vector, 0.0, b);

    // JACOBI
    jacobian_method(size, A, b, 0.0000000001);
    chebyshev_method(size, A, b, x, 1000, 0, 100000);

    free_vetor(x);
    free_vetor(A_1d);
    free_matrix(size, A);
    free(b);

    return 0;
}