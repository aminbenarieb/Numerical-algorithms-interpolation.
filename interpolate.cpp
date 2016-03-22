#include <stdlib.h>
#include <stdio.h>

#include "interpolate.h"
#include "inout.h"

int find_x(double* args, int size, double x){
	if ((args[0] <= args[size-1]) && x < args[0])
		return -1;

	int pos=0;
	for (int i=1; i<size; i++)
	{
		if ((args[i-1] <= x && x < args[i]) || //возрастающая таблица
			(args[i-1] >= x && x > args[i]))   //убывающая таблица
			break;
		pos++;
	}
	
	return pos;
}

double newton_polynom(double *args, double *funcs, int size, double x, int n, int *code){
    
    double **nodes = (double**) malloc(2*sizeof(double*)); //таблица x0,y0; x1,y1;...xn,yn
    for (int k=0; k<2; k++)
        nodes[k] = (double*) malloc((n+1)*sizeof(double));
    
    int pos = find_x(args, size, x); //позиция, относительно которой собирать x0,y0...xn,yn
    
#pragma region получаем массив узлов x0,y0
    int i=0; //в какой элемент узлов записываем
    int j=0; //из какого элемента таблицы
    //проходимся по всем элементам таблицы, в случае чего - разорвем цикл
    for (int k=0; ; k++)
    {
        //для позиции 5 из 7, вначале получим args[5], потом [6], [4], [7], [3], [2],...
        if (k%2==0)
            pos -= j;
        else
            pos += j;
        
        //проверяем на экстраполяцию
        if (pos >= size)
        {
            pos = size-1;
            *code = 1;
            continue;
        }
        else if (pos<0)
        {
            pos = 0;
            *code = 1;
            continue;
        }
        
        //наконец записываем число
        nodes[0][i] = args[pos];
        nodes[1][i] = funcs[pos];
        i++;
        j++;
        if (i == n+1)
            break;
    }
#pragma endregion
    
    //puts("\n[Nodes]:");
    //print_table(nodes,n+1,2);
    //puts("");
    
#pragma region получаем массив разделенных разностей
    /*разности хранятся в треугольном виде, как нормальный массив dd[строка][столбец]
     y01   y12  y23
     y012  y123
     y0123*/
    double **DivDif = (double**) malloc((n)*(n)*sizeof(double));
    for (int k=0; k<n+1; k++)
    {
        DivDif[k] = (double*) malloc((n)*sizeof(double));
        for (int j=0; j<n+1; j++)
            DivDif[k][j] = 0.0;
    }
    
    //инициализируем весь первый ряд yi,i+1
    //y(01) = (y0-y1)/(x0-x1)
    for (int i=0; i<n; i++)
    {
        DivDif[0][i] = (nodes[1][i]-nodes[1][i+1]) / (nodes[0][i]-nodes[0][i+1]);
    }
    //инициализируем остальные ряды через предыдущие
    for (int j=1; j<n; j++) //спускаемся по строкам, 01, 012, 0123
        for (int i=0; i<n-j; i++) //смещаемся по столбцам, 012, 123, 234
            DivDif[j][i] = (DivDif[j-1][i]-DivDif[j-1][i+1]) / (nodes[0][i]-nodes[0][i+j+1]);
#pragma endregion
    //PrintTable(DivDif, n, n);
    
    //наконец-то начинаем интерполировать полиномом лагранжа
    //Pn (x)=y0 + (x-x0)y01 + (x-x0)(x-x1)y012 + (x-x0)(x-x1)(x-x2)y0123
    double result=nodes[1][0]; //=y0
    
    for (int i=0; i<n; i++)
    {
        double mnozh = DivDif[i][0]; // =y012
        for (int j=0; j<=i; j++)
        {
            mnozh *= (x-nodes[0][j]);
        }
        result += mnozh;
    }
    
    return result;
}

double spline(double *args, double *funcs, int n, double x, int &code){
    
    
	int pos = find_x(args,n,x);
	if (pos < 0 || pos >= n)
	{
		code = 1;
		return 0.0;
	}
	pos += 1;

	//коэффициенты полинома
	double	*Ca = new double[n+2],
			*Cb = new double[n+2],
			*Cc = new double[n+2],
			*Cd = new double[n+2];
    
	double	*h = new double[n+1]; //шаг между соседними х
    
	//параметры трёхдиагональной матрицы; Ai y_(i-1) - Bi yi + Di yi = -Fi
	double	*A = new double[n+2], // Ai = h_(i-1)
			*B = new double[n+2], // Bi = -2(h_(i-1) + hi)
			*D = new double[n+2], // Di = hi
			*F = new double[n+2]; // Fi = -3 * ( (yi-y_(i-1))/hi - (y_(i-1)-y_(i-2))/h_(i-1) )
    
    
	//yi === ci = alpha_(i+1) c_(i+1) + beta_(i+1)
	double *alpha = new double[n+2]; //прогоночный Кси
	double *beta = new double[n+2]; //прогоночный Эта

	//определяем шаги между точками в x#
	for (int i=1; i<=n; i++)
		h[i] = args[i] - args[i-1];

	//инициализируем стартовые значения коэффициентов
	for (int i=1; i<=n; i++)
		Ca[i] = funcs[i-1]; //ai = y_(i-1)
	Cc[1] = Cc[n+1] = 0; //на границах 2 производная ==0
	alpha[2] = beta[2] = 0; //c1=0 => 0 = alpha2 c2 + beta2

	//в прогонке считать мы начнем с третьей альфы
	for (int i=2; i<=n; i++)
	{
		A[i] = h[i-1];
		B[i] = -2* (h[i-1] + h[i]);
		D[i] = h[i];
		F[i] = -3* ((funcs[i]-funcs[i-1])/h[i] - (funcs[i-1]-funcs[i-2])/h[i-1]);
	}

	//прогонка, прямой ход: определяем альфы и беты
	for (int i=2; i<=n; i++)
	{
		alpha[i+1] = D[i] / (B[i] - A[i]*alpha[i]);
		beta[i+1] = (A[i]*beta[i]+F[i]) / (B[i] - A[i]*alpha[i]);
	}
	
	//прогонка, обратный ход: определяем с
	for (int i=n; i>1; i--)
		Cc[i] = alpha[i+1]*Cc[i+1] + beta[i+1]; // c_(n+1)==0 => cn = 0 + beta_(n+1)

	//определяем b и d
	for (int i=1; i<=n; i++)
	{
		Cb[i] = (funcs[i]-funcs[i-1])/h[i] - h[i]/3 * (2*Cc[i] + Cc[i+1]);
		Cd[i] = (Cc[i+1] - Cc[i])/(3*h[i]);
	}
	//F(x) = ai + bi(x-x_(i-1)) + ci(x-x_(i-1))^2 + di(x-x_(i-1))^3
	double res = Ca[pos] + Cb[pos]*(x-args[pos-1]) + Cc[pos]*(x-args[pos-1])*(x-args[pos-1]) + Cd[pos]*(x-args[pos-1])*(x-args[pos-1])*(x-args[pos-1]);
	delete [] Ca, Cb, Cc, Cd, h, A, B, D, F, alpha, beta;
	return res;
}

double multi(double** matr, int sizeX, int sizeY, double x, int nX, double y, int nY, int &code){
	code = 0;

    
    // Берем строку X и столбец Y для полинома Ньютона
	double *XCol = new double[sizeX];
	for (int i=1; i<=sizeX; i++)
		XCol[i-1] = matr[i][0];

	double *YRow = new double[sizeY];
	for (int i=1; i<=sizeY; i++)
		YRow[i-1] = matr[0][i];

	double *ZRow = new double[sizeY];
    double *ZShear = new double[sizeX];
	for (int i=1; i<=sizeX; i++)
	{
		for (int j=1; j<=sizeY; j++)
			ZRow[j-1] = matr[i][j];
        
        ZShear[i-1] =  newton_polynom(YRow, ZRow, sizeY, y, nY, &code); //spline(YRow,ZRow,sizeY,y,code);
	}
    

	//срез собран - интерполируем его по Х
    
	double res = newton_polynom(XCol, ZShear, sizeX, x, nX, &code); //spline(XCol,ZShear,sizeX,x,code);
	delete [] XCol, YRow, ZRow, ZShear;

	return res;
}
