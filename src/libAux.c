#include<time.h>
#include<sys/time.h>
#include<stdio.h>
#include<stdlib.h>
#include<likwid.h>
#include<likwid-marker.h>
#include"libAux.h"


/* 
    Função que retorna tempo em milisegundos desde EPOCH
*/
double timestamp()
{
  struct timespec tp;

  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);

  return ( (double) tp.tv_sec*1.0e3 + (double) tp.tv_nsec*1.0e-6 );
}


/*
    Função que aloca um matriz
    Parametros:
        n: Ordem da matriz
*/
double **createMatrix(int n)
{
    double **M;

    /* Aloca as colunas */
    M = (double**) malloc(n*sizeof(double*));
    /* Verificar se memória foi alocada */
    if(M == NULL)
        exit(1);

    /* Aloca as linhas */
    for(int i = 0; i < n; i++)
    {
        M[i] = (double*) malloc(n*sizeof(double));
        /* Verifica se memória foi alocada */
        if(M[i] == NULL)
        {
            /* Libera posição já alocadas*/
            for(int j = 0; j < i; j++)
                free(M[j]);
            /* Libera primeira alocação */
            free(M);
            exit(1);
        }
    }
    return M;
}

/*
    Função que libera a memória alocada para uma matriz
    Parametros:
        M: Matriz
        n: Ordem da matriz
*/
double **destroyMatrix(double **M, int n)
{
    /* Libera as linhas */
    for(int i = 0; i < n; i++)
        free(M[i]);
    /* Libera as colunas */
    free(M);
    return NULL;
}

/*
    Função que copia uma matriz para outra
    Paramentros:
        original: Matriz original
        copy: Matriz cópia
        n: Ordem da matriz
*/
void copyMatrix(double **original, double **copy, unsigned int n)
{
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            copy[i][j] = original[i][j];
}

/*
    Função que imprime os elementos de uma matriz
    Paramentros:
        n: Ordem da matriz
        A: Matriz
*/
void printMatrix(int n, double **A)
{
    /* Imprime a matriz */
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            printf("%lf\t ", A[i][j]);
        printf("\n");
    }
}

/*
    Função que lê os elementos de uma matriz e do vetor de termos independentes
    Parametros: 
        n: Ordem da matriz
        A: Matriz
*/
void readLinearSystem(int n, double **A, double **Copy, double *b, double *cb)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(scanf("%lf", &A[i][j]) != 1)
            {
                printf("Erro ao ler a matriz!\n");
                exit(1);
            }
            Copy[i][j] = A[i][j];
        }
        if(scanf("%lf", &b[i]) != 1)
        {
            printf("Erro ao ler o vetor de termos independentes!\n");
            exit(1);
        }
        cb[i] = b[i];
    }
}

/*
    Função que imprime os elementos de uma matriz e do vetor de termos independentes
    Parametros: 
        n: Ordem da matriz
        A: Matriz
*/
void printLinearSystem(int n, double **A, double *b)
{
    /* Imprime a matriz */
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            printf("%lf\t ", A[i][j]);
        printf(" = ");
        printf("\t%lf\n", b[i]);
    }
}

/*
    Função que copia um vetor
    Parametros:
        original: Vetor original
        n: Ordem do vetor
*/
void copyVector(double *original, double*copy, unsigned int n)
{
    for(int i = 0; i < n; i++)
        copy[i] = original[i];
}

/*
    Função que imprime os elementos de um vetor
    Parametros: 
        n: Ordem da matriz
        x: Vetor a ser impresso
*/
void printVector(int n, double *x)
{
    /* Imprime a matriz */
    // for(int i = 0; i < n; i++)
        // printf("x[%d] = %lf\n", i, x[i]);
    printf("[ ");
    for(int i = 0; i < n; i++)
        printf("%lf ", x[i]);  
    printf("]\n");
}

/*
    Função que faz o cálculo do vetor residual
    Parametros:
        A: Matriz de coeficientes
        b: Vetor de valores independentes
        x: Vetor solução
        r: Vetor residual
        n: Ordem da matriz
*/
void residualVector(double **A, double *b, double *x, double *r, unsigned int n)
{
    for(int i = 0; i < n; i++)
    {
        double soma = 0.00; 
        for(int j = 0; j < n; j++)
            soma = soma + A[i][j]*x[j];
        soma = soma - b[i];
        r[i] = soma;
    }
}


/*
    Função que resolve sistemas triangulares superiores
    Parametros:
        A: Matriz
        b: Vetor de termos independentes
        x: Vetor solução
        n: Ordem da matriz
*/
void retroSubstitution(double **A, double *b, double *x, unsigned int n)
{
    /* Começa na última linha */
    LIKWID_MARKER_START("retroSubstitution");      
    for(int i = n-1; i >= 0; --i)
    {
        /* Atribui o valor da solução*/
        x[i] = b[i];

        /* Vai para a posição de um indíce acima e subtrai
         * os valores já calculados */
        for(int j = i+1; j < n; ++j)
            x[i] = x[i] - A[i][j]*x[j];

        /* Divide pelo coeficiente*/
        if( A[i][i] != 0.0)
            x[i] = x[i] / A[i][i];
        else
            x[i] = 0.0;
    }
    LIKWID_MARKER_STOP("retroSubstitution");
}

/*
    Função que faz a troca de duas linhas de um sistema linear
    Parametros:
        A: Matriz
        b: Vetor de termos independentes
        l1: Linha 1
        l2: Linha 2
*/
void swapLines(double **A, double *b, unsigned int l1, unsigned int l2)
{
    /* Troca as linhas da matriz */
    double *line = A[l1];
    A[l1] = A[l2];
    A[l2] = line;

    /* Troca os termos independentes */
    double ind = b[l1];
    b[l1] = b[l2];
    b[l2] = ind;
}

/*
    Função que encontrará o maior elemento do sistema linear
    Parametros:
        A: Matriz
        l: Linha atual
*/
unsigned int findMax(double **A, unsigned int l, unsigned int n)
{
    /* Considera o primeiro termo como o maior */
    double max = A[l][l];

    /* O pivo está na primeira linha */
    unsigned int lPivo = l;

    /* Procura por um pivô maior */
    for(int i = l+1; i < n ; i++)
    {
        if( A[i][l] > max )
        {
            max = A[i][l];
            lPivo = i;
        }
    }

    /* Retorna em qual linha está o pivo */
    return lPivo;
}

/*
    Função que aplica a Eliminação-Gaussiana com pivoteamento parcial
    Parametros:
        A: Matriz
        b: Vetor de termos independentes
        n: Ordem da matriz
*/ 
void classicEliminationWithPivot(double **A, double *b, unsigned int n)
{
    /* Inicia na primeira linha e primeira coluna */
    for(int i = 0; i < n ; ++i)
    {
         /* Pivoteamento parcial */
        unsigned int iPivot = findMax(A,i,n);
        if( i != iPivot )
            swapLines(A,b,i,iPivot);

        if( A[i][i] == 0.0 )
        {   
            int j = i+1;
            while( A[j][i] == 0.0 && j < n-1 )
                j++;
            if( j == n-1 )
            {
                printf("Coluna com todos coeficientes igual a zero\n");
                exit(1);
            }
            swapLines(A,b,i,j);
        }
        
        LIKWID_MARKER_START("classicEliminationWithPivot");
        /* Vai para a próxima linha */
        for(int k = i+1; k < n; k++)
        {
            /* Calcula o multiplicador */
            double m = A[k][i] / A[i][i];

            /* Para garantir que o elemento abaixo do pivo será zero */
            A[k][i] = 0.0;

            /* Subtrai de cada elemento da equação a multiplicação
             * desde elemento com o multiplicador calculado */
            for(int j = i+1; j < n; ++j)
                A[k][j] = A[k][j] - (m * A[i][j]);
            b[k] = b[k] - m*b[i];
        }
        LIKWID_MARKER_STOP("classicEliminationWithPivot");
    }
}

/*
    Função que aplica a Eliminação-Gaussiana sem pivoteamento parcial
    Parametros:
        A: Matriz
        b: Vetor de termos independentes
        n: Ordem da matriz
*/
void classicEliminationWithoutMult(double **A, double *b, unsigned int n)
{
    /* Inicia na primeira linha e primeira coluna */
    for(int i = 0; i < n ; ++i)
    {
        LIKWID_MARKER_START("classicEliminationWithPivotWithoutMult");
        /* Pivoteamento parcial */

        unsigned int iPivot = findMax(A,i,n);
        if( i != iPivot )
            swapLines(A,b,i,iPivot);

        if( A[i][i] == 0.0 )
        {   
            int j = i+1;
            while( A[j][i] == 0.0 && j < n-1 )
                j++;
            if( j == n-1 )
            {
                printf("Coluna com todos coeficientes igual a zero\n");
                exit(1);
            }
            swapLines(A,b,i,j);
        }
        
        /* Vai para a próxima linha */
        for(int k = i+1; k < n; ++k)
        {
            for(int j = i+1; j < n; ++j)
                A[k][j] = A[k][j]*A[i][i] - A[i][j]*A[k][i];
                
            b[k] = b[k]*A[i][i] - b[i]*A[k][i];
            A[k][i] = 0.0;
        }
        LIKWID_MARKER_STOP("classicEliminationWithPivotWithoutMult");
    }
}

/*  > Para cada linha  do pivô, dividir todos os  coeficientes pelo pivô  (que fica com o valor 1)
    > Proceder com  a eliminação,  zerando a coluna  do pivô,  sem fazer pivoteamento.
    > Completada   a  triangularização,   calcular  as   incógnitas  por retro-substituição.
*/
/*
    Função que aplica a forma alternativa de eliminação
    Parametros:
        A: Matriz
        b: Vetor de termos independentes
        n: Ordem da matriz
*/
void alternativeFormOfElimination(double **A, double *b, unsigned int n)
{
    LIKWID_MARKER_START("alternativeFormOfElimination");
    /* Inicia na primeira linha e primeira coluna */
    for(int i = 0; i < n ; ++i)
    {    
        double pivo = A[i][i];
        if( pivo == 0.0 )
        {
            int j = i+1;
            while( A[j][i] == 0.0 && j < n-1 )
                j++;
            if( j == n-1 )
            {
                printf("Coluna com todos coeficientes igual a zero\n");
                exit(1);
            }
            swapLines(A,b,i,j);
            pivo = A[i][i];
        }

        /* Divide cada elemento da linha do pivô, pelo pivô */
        for(int l = 0; l < n; ++l)
            A[i][l] = A[i][l] / pivo;
        b[i] = b[i] / pivo;   
        
        /* Vai para a próxima linha */
        for(int k = i+1; k < n; ++k)
        {
            /* Calcula o multiplicador */
            double m = A[k][i] / A[i][i];
            
            /* Para garantir que o elemento será zero */
            A[k][i] = 0.0;

            /* Subtrai de cada elemento da equação a multiplicação
             * desde elemento com o multiplicador calculado */
            for(int j = i+1; j < n; ++j)
                A[k][j] = A[k][j] - m*A[i][j];
            b[k] = b[k] - m*b[i];
        }
    } 
    LIKWID_MARKER_STOP("alternativeFormOfElimination");
}

