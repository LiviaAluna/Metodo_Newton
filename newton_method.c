#include<stdio.h>
#include<math.h>
#define H 0.000001
#define N 2
#define M 4
#define O 1

double funcao(double x, double y){
    double f; 
    f=x*x+y*y;
    return f;
}

double derX(double x, double y){
    double dx;
    dx=(funcao(x+H,y)-funcao(x,y))/H;
    return dx;
}
double derY(double x, double y){
    double dy;
    dy=(funcao(x,y+H)-funcao(x,y))/H;
    return dy;
}
double derX2(double x, double y){
    double dxx;
    dxx=(funcao(x+H,y)-2*funcao(x,y)+funcao(x-H,y))/(H*H);
    return dxx;
}
double derY2(double x, double y){
    double dyy;
    dyy=(funcao(x,y+H)-2*funcao(x,y)+funcao(x,y-H))/(H*H);
    return dyy;
}
double derXY(double x, double y){
    double dxy;
    dxy=(funcao(x+H,y+H)-funcao(x+H,y-H)-funcao(x-H,y+H)+funcao(x-H,y-H))/(4*(H*H));
    return dxy;
}
void multMatrix(double A[][N], double B[][N], double C[][N])
{
    // Percorre as linhas de A
    for (int i = 0; i < N; i++){
        // Percorre as linhas de B
        for (int j = 0; j < N; j++){
            C[i][j] = 0.0;
            // Percorre a linha de A e a coluna de B
            for(int k = 0; k < N; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void printMatrix(double matrix[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.1f\t", matrix[i][j]);
        }
        printf("\n");
    }
}

void main(){
    double x = 2, y= 2;
//Matriz Jacobiana

    double jac[N][O];
    jac[0][0]=derX(x,y);
    jac[1][0]=derY(x,y);
    printf("\nmatriz jacobiana:\n%.1f, %.1f\n", jac[0][0],jac[1][0]);

//Norma Vetorial

    double soma = 0.0;
    for(int i = 0; i < N; i++){    
    soma +=jac[i][0]*jac[i][0];
    
    }

    double norma = sqrt(soma);

    printf("norma vetorial=%.1f\n", norma);

//Matriz Hessiana

    double hes[N][N];
    hes[0][0]=derX2(x,y);
    hes[0][1]=derXY(x,y);
    hes[1][0]=derXY(x,y);
    hes[1][1]=derY2(x,y); 
    printf("\nmatriz hessiana:\n%.1f  %.1f\n%.1f  %.1f\n", hes[0][0], hes[0][1], hes[1][0], hes[1][1]);

//Matriz inversa da Hessiana
    double d,adjHes[2][2];
    int i,j;
    double invHes[2][2];

    /* Calcular o determinante da matriz A */
    d = (hes[0][0]*hes[1][1])-(hes[0][1]*hes[1][0]);
    if(d==0)
    {
        printf("Determinante 0");
        return 0;
    }
     /* Encontrar adjHes */
    adjHes[0][0]=hes[1][1];
    adjHes[1][1]=hes[0][0];
    adjHes[0][1]=-hes[0][1];
    adjHes[1][0]=-hes[1][0];

    /* Encontrar matriz inversa*/
    printf("Matriz inversa da hessiana: \n");
    for(i=0;i<2;i++)
    {
        for(    j=0;j<2;j++)
        {
            invHes[i][j]=(adjHes[i][j])/d;
        }
    }
    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            printf("%.1f ",invHes[i][j]);
        }
        printf("\n");
    }

//multiplicação das matrizes
   double resultMult[N][O];
   multMatrix(invHes,jac,resultMult);
   printf("Multiplicacao entre matriz hessiana inversa e matriz jacobiana:\n");
   for(i=0;i<2;i++)
    {
        for(j=0;j<1;j++)
        {
            printf("%.1f ",resultMult[i][j]);
        }
        printf("\n");
    }

}