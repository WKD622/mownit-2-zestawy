#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_interp.h>
#define N 10

int main(void)
{
  float f = 1.0 / 3.0;
  double d = 1.0 / 3.0;

  double fd = f; /* promote from float to double */

  // printf(" f="); gsl_ieee_printf_float(&f);
  // printf("\n");

  // printf("fd="); gsl_ieee_printf_double(&fd);
  // printf("\n");

  int i;
  for (i = 0; i < N; i++)
  {
    printf("%d\n d=", i + 1);
    gsl_ieee_printf_double(&d);
    printf("\n");
    printf(" d= %.50f\n", d);
    d = d / 4.5;
  }

  return 0;
}
