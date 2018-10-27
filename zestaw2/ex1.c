#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>
#define N 1000

double *generate_matrix_array(int size1, int size2)
{
  int i;
  double *array = malloc(size1 * size2 * sizeof(double));
  if (array == NULL)
    return NULL;
  for (i = 0; i < size1 * size2; i++)
  {
    array[i] = (double)rand() / RAND_MAX * 100.0;
  }
  return array;
}

double *generate_vector_array(int size)
{
  int i;
  double *array = malloc(size * sizeof(double));
  if (array == NULL)
    return NULL;
  for (i = 0; i < size; i++)
  {
    array[i] = (double)rand() / RAND_MAX * 100.0;
  }
  return array;
}

void print_matrix_array(double *matrix, int size1, int size2)
{
  int i;
  for (i = 0; i < size1 * size2; i++)
    printf("%f ", matrix[i]);
  printf("\n\n");
}

void print_vector_array(double *vector, int size)
{
  int i;
  for (i = 0; i < size; i++)
    printf("%f ", vector[i]);
  printf("\n\n");
}

void print_equation_1(double *matrix, double *vector, int size, char *title)
{
  printf("| %s\n", title);
  int i, j = 0;
  for (i = 0; i < size; i++)
  {
    if (i / size == 0)
      printf("|");
    for (j = size * i; j < size * (i + 1); j++)
    {
      printf("%10f ", matrix[j]);
    }
    if (i == size / 2)
      printf(" * ");
    else
      printf("   ");
    printf("  x%d  ", i + 1);
    if (i == size / 2)
      printf(" = ");
    else
      printf("   ");
    printf("%10f", vector[i]);
    printf("\n");
  }
  printf("\n");
}

void print_equation_2(double *matrix, gsl_vector *vector, int size, char *title)
{
  printf("| %s\n", title);
  int i, j = 0;
  for (i = 0; i < size; i++)
  {
    if (i / size == 0)
      printf("|");
    for (j = size * i; j < size * (i + 1); j++)
    {
      printf("%10f ", matrix[j]);
    }
    if (i == size / 2)
      printf(" * ");
    else
      printf("   ");
    printf("  %10f  ", vector->data[i]);
    if (i == size / 2)
      printf(" = ");
    else
      printf("   ");
    printf("  ?");
    printf("\n");
  }
  printf("\n");
}

void print_vector(gsl_vector *result, int size, char *title)
{
  int i;
  printf("| %s\n", title);
  for (i = 0; i < size; i++)
  {
    printf("| x%d", i + 1);
    if (i == size / 2)
      printf("  = ");
    else
      printf("    ");
    printf("%20.15f\n", result->data[i]);
  }
  printf("\n");
}

void print_vector_q(gsl_vector *result, int size, char *title)
{
  int i;
  printf("| %s\n", title);
  for (i = 0; i < size; i++)
  {
    printf("| ?");
    if (i == size / 2)
      printf("  = ");
    else
      printf("    ");
    printf("%20.15f\n", result->data[i]);
  }
  printf("\n");
}

void print_diffrence(gsl_vector *result, double *B, int size)
{
  int i;
  printf("| DIFFRENCE:\n");
  for (i = 0; i < size; i++)
  {
    printf("| %24.20f \n", result->data[i] - B[i]);
  }
  printf("\n");
}

void solve(int size, gsl_matrix_view m, gsl_vector_view b, gsl_vector *x, gsl_permutation *p)
{
}

double czas(struct rusage *ru0, struct rusage *ru1, FILE *result)
{

  double utime = 0, stime = 0, ttime = 0;

  /* Obliczenie czasow. Aby mikrosekundy traktowac jako czesci sekund musimy je wymnozyc przez 10^-6*/
  utime = (double)ru1->ru_utime.tv_sec + 1.e-6 * (double)ru1->ru_utime.tv_usec - ru0->ru_utime.tv_sec - 1.e-6 * (double)ru0->ru_utime.tv_usec;
  stime = (double)ru1->ru_stime.tv_sec + 1.e-6 * (double)ru1->ru_stime.tv_usec - ru0->ru_stime.tv_sec - 1.e-6 * (double)ru0->ru_stime.tv_usec;
  ttime = stime + utime;

  /*printf("user time: %3f\n", utime);
	printf("system time: %3f\n", stime);
  printf("total time: %3f\n", ttime);*/
  return ttime;
}

void test()
{
  int size;
  struct rusage t0, t1, t2;
  int s;
  FILE *test_result = fopen("out/test_result", "wr");
  for (size = 10; size < 1000; size++)
  {
    double *B = generate_vector_array(size);
    double *A = generate_matrix_array(size, size);

    gsl_matrix_view m = gsl_matrix_view_array(A, (size_t)size, (size_t)size);
    gsl_vector_view b = gsl_vector_view_array(B, (size_t)size);

    gsl_permutation *p = gsl_permutation_alloc((size_t)size);
    gsl_vector *x = gsl_vector_alloc((size_t)size);

    getrusage(RUSAGE_SELF, &t0);
    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    getrusage(RUSAGE_SELF, &t1);
    gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);
    getrusage(RUSAGE_SELF, &t2);

    system("clear");
    printf("| TESTING: \n| %d of 1000\n", size);
    double ttime1 = czas(&t0, &t1, test_result);
    double ttime2 = czas(&t1, &t2, test_result);
    fprintf(test_result, "%3f %3f\n", ttime1, ttime2);
  }
  fclose(test_result);
  system("clear");
  printf("| DONE:\n| Results are in out folder.\n");
  system("gnuplot --persist -e 'plot \"out/test_result\" u 1:2'");
  sleep(3);
  system("clear");
}

int main(int argc, char *argv[])
{
  srand(time(NULL));

  if (strcmp(argv[1], "-t") == 0)
    test();
  else
  {
    int size = atoi(argv[1]);
    int s;

    double *B = generate_vector_array(size);
    double *A = generate_matrix_array(size, size);
    double *A_copy = malloc(size * size * sizeof(double));
    double *B_copy = malloc(size * sizeof(double));

    memcpy(A_copy, A, size * size * sizeof(double));
    memcpy(B_copy, B, size * sizeof(double));
    print_equation_1(A, B, size, "START:");

    gsl_matrix_view m = gsl_matrix_view_array(A, (size_t)size, (size_t)size);
    gsl_matrix_view m_copy = gsl_matrix_view_array(A_copy, (size_t)size, (size_t)size);
    gsl_vector_view b = gsl_vector_view_array(B, (size_t)size);
    gsl_vector_view b_copy = gsl_vector_view_array(B_copy, (size_t)size);

    gsl_permutation *p = gsl_permutation_alloc((size_t)size);
    gsl_vector *x = gsl_vector_alloc((size_t)size);
    gsl_vector *result = gsl_vector_alloc((size_t)size);

    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);

    print_vector(x, size, "RESULT:");

    print_equation_2(A_copy, x, size, "LET'S CHECK:");

    gsl_blas_dgemv(CblasNoTrans, 1.0, &m_copy.matrix, x, 0.0, result);

    print_vector_q(result, size, "RESULT:");
    print_diffrence(result, B_copy, size);

    gsl_permutation_free(p);
    gsl_vector_free(x);
    gsl_vector_free(result);
    free(A);
    free(B);
    free(A_copy);
  }
  return 0;
}
