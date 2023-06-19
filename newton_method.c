#include<stdio.h>
#include<math.h>
#define N 2
#define O 1

double hessiana[N][N],jacobiana[N][O],invH[N][N],adjHes[N][N],x, y,e=0.0001, h=0.00000001,d[N][O];

double funcao(double x, double y){
    double f; 
    f=x*x+y*y;
    return f;
}
double derX(double x, double y){
    double dx;
    dx=(funcao(x+h, y)-funcao(x-h,y))/(2*h);
    return dx;
}
double derY(double x, double y){
    double dy;
    dy=(funcao(x,y+h)-funcao(x,y-h))/(2*h);
    return dy;
}
double derX2(double x, double y){
    double dxx;
    dxx=(funcao(x+h,y)-2*funcao(x,y)+funcao(x-h,y))/(h*h);
    return dxx;
}
double derY2(double x, double y){
    double dyy;
    dyy=(funcao(x,y+h)-2*funcao(x,y)+funcao(x,y-h))/(h*h);
    return dyy;
}
double derXY(double x, double y){
    double dxy;
    dxy=(funcao(x+h,y+h)-funcao(x+h,y-h)-funcao(x-h,y+h)+funcao(x-h,y-h))/(4*(h*h));
    return dxy;
}
void invhess(){
    double det,adjHes[N][N];
    int i,j;

    /* Calcular o determinante da matriz A */
    det = (hessiana[0][0]*hessiana[1][1])-(hessiana[0][1]*hessiana[1][0]);
    if(det==0)
    {
        return 0;
    }
     /* Encontrar adjHes */
    adjHes[0][0]=hessiana[1][1];
    adjHes[1][1]=hessiana[0][0];
    adjHes[0][1]=-hessiana[0][1];
    adjHes[1][0]=-hessiana[1][0];

    /* Encontrar matriz inversa*/
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            invH[i][j]=(adjHes[i][j])/det;
        }
    }
}
void hess(){
    hessiana[0][0]=derX2(x,y);
    hessiana[0][1]=derXY(x,y);
    hessiana[1][0]=derXY(x,y);
    hessiana[1][1]=derY2(x,y);
}
void jac(){
    jacobiana[0][0]=derX(x,y);
    jacobiana[1][0]=derY(x,y);
}
double norma(double x[][O]){
    int i;
    double acum= 0;

    for(i=0; i<N; i++){
        acum+=x[i][O]*x[i][O];
    }
    acum=sqrt(acum);

    return acum;
}
void main(){
    x = 1.0; y= 1.0;
    int i, j, it=0, itMax=100;
    jac();
    //imprimir jacobiana(gradiente) na tela
    for(i=0;i<N;i++) for(j=0;j<O;j++) printf("jacobiana[%d][%d]=%f\n",i,j, jacobiana[i][j]);
    //imprimir hessiana na tela
    hess();
    for(i=0;i<N;i++) for(j=0; j<N;j++) printf("H[%d][%d]=%f\n",i,j, hessiana[i][j]);
    //Printar a matriz hessiana inversa
    invhess();
    for(i=0;i<N;i++) for(j=0; j<N;j++) printf("invH[%d][%d]=%f\n",i,j, invH[i][j]);
    //Encontrar o minimo
    while(norma(jacobiana)>e && it<itMax){
        d[0][0]=(-1)*(invH[0][0]*jacobiana[0][0]+invH[0][1]*jacobiana[1][0]);
        d[1][0]=(-1)*(invH[1][0]*jacobiana[0][0]+invH[1][1]*jacobiana[1][0]);
            
         x+=d[0][0];
         y+=d[1][0];

        jac();
        hess();
        invhess();

       printf("interacao = %d \t f(%f)(%f)=%f\n",it, x,y, funcao(x, y));
    it++;
    }

}
