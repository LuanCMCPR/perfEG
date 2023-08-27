#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"libAux.h"

int main()
{
    FILE *fp;
    clock_t s, e;
    int n;
    double **A;
    double *b, *x, *r;
    
    // printf("Digite a ordem da matriz: ");
    scanf("%d", &n);

    /* Aloca a matriz */
    A = createMatrix(n);
    /* Aloca o vetor b */
    b = (double*) malloc(n*sizeof(double));
    /* Aloca o vetor x */
    x = (double*) malloc(n*sizeof(double));
    /* Aloca o vetor r */
    r = (double*) malloc(n*sizeof(double));


    // printf("Digite os elementos da matriz: \n");   
    /* Lê a matriz */
    readLinearSystem(n, A, b);

    printf("\nSistema linear Digitado: \n");
    /* Imprime o sistema linear */
    printLinearSystem(n, A, b);


    /***** Forma clássica com pivo *****/
    // printf("\nForma clássica com pivo: \n");
    // classicEliminationWithPivot(A, b, n);
    // /* Resolve o sistema triangular superior */
    // retroSubstituion(A, b, x, n);
    // /* Imprime o vetor solução */
    // printf("Solução: \n");
    // printVector(n, x);
    // /* Calcula o vetor residual */
    // printf("Resíduo: \n");
    // residualVector(A, b, x, r, n);
    // /* Imprime o vetor residual */
    // printVector(n, r);


    
    /***** Forma clássica com pivo e sem multiplicador *****/
    // printf("\nForma clássica com pivo e sem multiplicador: \n");
    // classicEliminationWithoutMult(A, b, n);
    // /* Resolve o sistema triangular superior */
    // retroSubstituion(A, b, x, n);
    // /* Imprime o vetor solução */
    // printf("Solução: \n");
    // printVector(n, x);
    // /* Calcula o vetor residual */
    // printf("Resíduo: \n");
    // residualVector(A, b, x, r, n);
    // /* Imprime o vetor residual */
    // printVector(n, r);



    /***** Forma alternativa *****/
    printf("\nForma alternativa de eliminação: \n");
    alternativeFormOfElimination(A, b, n);
    /* Resolve o sistema triangular superior */
    retroSubstituion(A, b, x, n);
    /* Imprime o vetor solução */
    printf("Solução: \n");
    printVector(n, x);
    /* Calcula o vetor residual */
    printf("Resíduo: \n");
    residualVector(A, b, x, r, n);
    /* Imprime o vetor residual */
    printVector(n, r);

    free(x);
    free(b);
    free(r);

    return 0;

}