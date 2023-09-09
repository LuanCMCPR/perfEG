

/* Struct para um sistema linear */
typedef struct sistemaLinear
{
  int n; /* Order of the linear system */
  double **A; /* Matrix with coeficients */
  double *b; /* Independent Terms*/
} SL_t;

/* 
    Função que retorna tempo em milisegundos desde EPOCH
*/
double timestamp();

/* 
    Função que aloca um matriz
    Parametros:
        n: Ordem da matriz
 */
double **createMatrix(int n);

/*
    Função que libera a memória alocada para uma matriz
    Parametros:
        M: Matriz
        n: Ordem da matriz
*/
double **destroyMatrix(double **M, int n);

/*
    Função que copia uma matriz para outra
    Paramentros:
        original: Matriz original
        copy: Matriz cópia
        n: Ordem da matriz
*/
void copyMatrix(double **original, double **copy, unsigned int n);

/*
    Função que imprime os elementos de uma matriz
    Paramentros:
        n: Ordem da matriz
        A: Matriz
*/
void printMatrix(int n, double **A);

/*
    Função que lê os elementos de uma matriz e do vetor de termos independentes
    Parametros: 
        n: Ordem da matriz
        A: Matriz
*/
void readLinearSystem(int n, double **A, double **Copy, double *b, double *cb);

/*
    Função que imprime os elementos de uma matriz e do vetor de termos independentes
    Parametros: 
        n: Ordem da matriz
        A: Matriz
*/
void printLinearSystem(int n, double **A, double *b);

/*
    Função que imprime os elementos de um vetor
    Parametros: 
        n: Ordem da matriz
        x: Vetor
*/
void printVector(int n, double *x);

/*
    Função que copia um vetor
    Parametros:
        original: Vetor original
        copy: Vetor cópia
        n: Ordem do vetor
*/
void copyVector(double *original, double*copy, unsigned int n);

/*
    Função que faz o cálculo do vetor residual
    Parametros:
        A: Matriz de coeficientes
        b: Vetor de valores independentes
        x: Vetor solução
        r: Vetor residual
        n: Ordem da matriz
*/
void residualVector(double **A, double *b, double *x, double *r, unsigned int n);


/*
    Função que resolve sistemas triangulares superiores
    Parametros:
        A: Matriz
        b: Vetor de termos independentes
        x: Vetor solução
        n: Ordem da matriz
*/
void retroSubstitution(double **A, double *b, double *x, unsigned int n);

/*
    Função que aplica a Eliminação-Gaussiana com pivoteamento parcial
    Parametros:
        A: Matriz
        b: Vetor de termos independentes
        n: Ordem da matriz
*/ 
void classicEliminationWithPivot(double **A, double *b, unsigned int n);

/*
    Função que aplica a Eliminação-Gaussiana sem pivoteamento parcial
    Parametros:
        A: Matriz
        b: Vetor de termos independentes
        n: Ordem da matriz
*/
void classicEliminationWithoutMult(double **A, double *b, unsigned int n);

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
void alternativeFormOfElimination(double **A, double *b, unsigned int n);
