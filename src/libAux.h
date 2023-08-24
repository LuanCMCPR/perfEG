

/* Struct para um sistema linear */
typedef struct sistemaLinear
{
  int n; /* Ordem of the linear system */
  double **A; /* Matrix with coeficients */
  double *b; /* Independent Terms*/
  double *x; /* Solution for de linear system */

} SL_t;

/* 
Função que aloca um matriz
Parametros:
    n: Ordem da matriz
 */
double **createMatrix(int n);

/*
Função que lê os elementos de uma matriz
Parametros: 
    n: Ordem da matriz
    A: Matriz
*/
void readMatrix(int n, double **A);


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
void readLinearSystem(int n, double **A, double *b);

/*
Função que imprime os elementos de uma matriz e do vetor de termos independentes
Parametros: 
    n: Ordem da matriz
    A: Matriz
*/
void printLinearSystem(int n, double **A, double *b);

/* Resolver o sistema linear */
void retroSubstituion(double **A, double *b, double *x, unsigned int n);

/* Calcula a forma clássica com o pivo */
void classicEliminationWithPivot();

/* Calcula a forma clássica sem o pivo */
void classicEliminationWithoutPivot();

/*  > Para cada linha  do pivô, dividir todos os  coeficientes pelo pivô  (que fica com o valor 1)
    > Proceder com  a eliminação,  zerando a coluna  do pivô,  sem fazer pivoteamento.
    > Completada   a  triangularização,   calcular  as   incógnitas  por retro-substituição.
*/
void alternativeFormOfElimination();
