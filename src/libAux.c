#include<stdio.h>
#include<stdlib.h>
#include"libAux.h"

/*
Função que aloca um matriz
Parametros:
    n: Ordem da matriz
 */
double **createMatrix(int n)
{
    double **M;

    /* Alocas as colunas */
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
Função para liberar memória alocada para um sistema linear
Parametros:
    SL: Sistema linear
*/
double destroyLinearSystem(SL_t *SL)
{
    /* Libera as linhas */
    for(int i = 0; i < SL->n; i++)
        free(SL->A[i]);
    /* Libera as colunas */
    free(SL->A);
    /* Libera o vetor b */
    free(SL->b);
    /* Libera o vetor x */
    free(SL->x);
    /* Libera o sistema linear */
    free(SL);
}

/*
Função que lê os elementos de uma matriz
Parametros: 
    n: Ordem da matriz
    A: Matriz
*/
void readMatrix(int n, double **A)
{

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            scanf("%lf", &A[i][j]);
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
void readLinearSystem(int n, double **A, double *b)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            scanf("%lf", &A[i][j]);
        scanf("%lf", &b[i]);
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
Função que imprime os elementos de um vetor
Parametros: 
    n: Ordem da matriz
    x: Vetor
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


void residualVector(double **A, double *b, double *x, double *r, unsigned int n)
{
    /* Calcula o vetor residual */
    // for(int i = 0; i < n; i++)
    // {
    //     r[i] = b[i];
    //     for(int j = 0; j < n; j++)
    //         r[i] = r[i] - A[i][j]*x[j];
    // }
    double soma = 0.00;
    for(int i = 0; i < n; i++)
    {
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
void retroSubstituion(double **A, double *b, double *x, unsigned int n)
{
    /* Começa na última linha */
    for(int i = n-1; i >= 0; i--)
    {
        /* Atribui o valor da solução*/
        x[i] = b[i];

        /* Vai para a posição de um indíce acima e subtrai
         * os valores já calculados */
        for(int j = i+1; j < n; ++j)
            x[i] = x[i] - A[i][j]*x[j];

        /* Divide pelo coeficiente*/
        x[i] = x[i] / A[i][i];
    }
}

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
    for(int l = 0; l < n ; l++)
    {
         /* Pivoteamento parcial */
        unsigned int lPivot = findMax(A,l,n);
        if( l != lPivot )
            swapLines(A,b,l,lPivot);
        
        /* Vai para a próxima linha */
        for(int i = l+1; i < n; i++)
        {
            /* Calcula o multiplicador */
            double m = A[i][l] / A[l][l];

            /* Para garantir que o elemento será zero */
            A[i][l] = 0.0;

            /* Subtrai de cada elemento da equação a multiplicação
             * desde elemento com o multiplicador calculado */
            for(int j = l+1; j < n; j++)
                A[i][j] = A[i][j] - m*A[l][j];
            b[i] = b[i] - m*b[l];
        }
    }
}

/*
Função que aplica a Eliminação-Gaussiana sem pivoteamento parcial
Parametros:
    A: Matriz
    b: Vetor de termos independentes
    n: Ordem da matriz
*/
void classicEliminationWithoutPivot(double **A, double *b, unsigned int n)
{
    /* Inicia na primeira linha e primeira coluna */
    for(int l = 0; l <n ; l++)
    {    
        /* Vai para a próxima linha */
        for(int i = l+1; i < n; i++)
        {
            /* Calcula o multiplicador */
            double m = A[i][l] / A[l][l];
            
            /* Para garantir que o elemento será zero */
            A[i][l] = 0.0;

            /* Subtrai de cada elemento da equação a multiplicação
             * desde elemento com o multiplicador calculado */
            for(int j = l+1; j < n; j++)
                A[i][j] = A[i][j] - m*A[l][j];
            b[i] = b[i] - m*b[l];
        }
    }

    /* Resolve o sistema */
}

/*  > Para cada linha  do pivô, dividir todos os  coeficientes pelo pivô  (que fica com o valor 1)
    > Proceder com  a eliminação,  zerando a coluna  do pivô,  sem fazer pivoteamento.
    > Completada   a  triangularização,   calcular  as   incógnitas  por retro-substituição.
*/
void alternativeFormOfElimination(double **A, double *b, unsigned int n)
{
    /* Inicia na primeira linha e primeira coluna */
    for(int l = 0; l <n ; l++)
    {    
        double pivo = A[l][l];

        /* Divide cada elemento da linha do pivô, pelo pivô */
        for(int i = 0; i < n; i++)
            A[l][i] = A[l][i] / pivo;
        b[l] = b[l] / pivo;   
        
        /* Vai para a próxima linha */
        for(int i = l+1; i < n; i++)
        {
            /* Calcula o multiplicador */
            double m = A[i][l] / A[l][l];
            
            /* Para garantir que o elemento será zero */
            A[i][l] = 0.0;

            /* Subtrai de cada elemento da equação a multiplicação
             * desde elemento com o multiplicador calculado */
            for(int j = l+1; j < n; j++)
                A[i][j] = A[i][j] - m*A[l][j];
            b[i] = b[i] - m*b[l];
        }
    } 
}

