#include<stdio.h>
#include<stdlib.h>
#include"libAux.h"

int main()
{
    int n;
    double **A;
    double *b, *x;
    
    printf("Digite a ordem da matriz: ");
    scanf("%d", &n);

    /* Aloca a matriz */
    A = createMatrix(n);
    /* Aloca o vetor b */
    b = (double*) malloc(n*sizeof(double));
    /* Aloca o vetor x */
    x = (double*) malloc(n*sizeof(double));

    printf("Digite os elementos da matriz: \n");   
    /* Lê a matriz */
    readLinearSystem(n, A, b);

    /* Imprime o sistema linear */
    printLinearSystem(n, A, b);

    retroSubstituion(A, b, x, n);

    printf("Solução do sistema: \n");
    for(int i = 0; i < n; i++)
        printf("x[%d] = %lf\n", i, x[i]);
    return 0;

}