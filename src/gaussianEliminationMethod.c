#include<time.h>
#include<sys/time.h>
#include<stdio.h>
#include<stdlib.h>
#include<likwid.h>
#include<likwid-marker.h>
#include"libAux.h"
#define MAX_METHODS 3


int main()
{
    int n;
    double **A, **CA, *b, *cb, *x, *r, diff;

    /* Lê a ordem da matriz */
    if(scanf("%d", &n) != 1)
    {
        printf("Erro ao ler a ordem da matriz!\n");
        exit(1);
    }

    /* Cria a matriz A e aloca os vetores b, r, x */    
    A = createMatrix(n);
    CA = createMatrix(n);

    b = (double*) malloc(n*sizeof(double));
    cb = (double*) malloc(n*sizeof(double));
    x = (double*) malloc(n*sizeof(double));
    r = (double*) malloc(n*sizeof(double));

    /* Lê a matriz original */
    readLinearSystem(n, A, CA, b, cb);

    LIKWID_MARKER_INIT;
    for (int i = 1; i <= MAX_METHODS; i++) 
    {
        diff = timestamp();
        switch (i) 
        {
            case 1:
                printf("\nFORMA CLÁSSICA COM PIVO\n");
                classicEliminationWithPivot(A, b, n);
                break;
            case 2:
                printf("\nFORMA CLÁSSICA COM PIVO E SEM MULTIPLICADOR\n");
                classicEliminationWithoutMult(A, b, n);
                break;  
            case 3:
                printf("\nFORMA ALTERNATIVA DE ELIMINAÇÃO\n");
                alternativeFormOfElimination(A, b, n);
                break;
        }

        /* Resolve o sistema triangular superior */
        retroSubstitution(A, b, x, n);

        /* Marca o fim da temporização, calcula a diferença e imprime o tempo de execução. */
        diff = timestamp() - diff;
        printf("Tempo de execução: %lf ms\n", diff);

        /* Imprime o vetor solução */
        printf("Solução: \n");
        printVector(n, x);

        /* Calcula o vetor residual */
        printf("Resíduo: \n");
        residualVector(A, b, x, r, n);
        printVector(n, r);

        /* Copia o sistema linear novamente */
        copyMatrix(CA, A, n);
        copyVector(cb, b, n);

        /* Zera vetores solução e residuo */
        for(int i = 0; i < n; i++)
        {
            x[i] = 0.00;
            r[i] = 0.00;
        }

    }
    
    /* Liberar memória alocada */
    free(x);
    free(cb);
    free(b);
    free(r);
    destroyMatrix(CA, n);
    destroyMatrix(A, n);

    LIKWID_MARKER_CLOSE;

    return 0;

}