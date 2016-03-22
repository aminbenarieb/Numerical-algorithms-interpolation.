//
//  main.c
//  lab1
//
//  Created by Amin Benarieb on 17/03/16.
//  Copyright © 2016 Amin Benarieb. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "inout.h"
#include "interpolate.h"

#define EULER 2.7182818284590452353602874713527
#define COLUMNS 2
#define SIZE 6
#define SPLINE_SIZE 10

double f(double x){
    return x*x+2;
}

void fill_table(double** &table, double left, double right, const int n){
    
    double step = ((right-left) / n);
    for (int i=0; i<n; i++)
    {
        table[0][i] = left + step*i;
        table[1][i] = f(table[0][i]);
    }
}

int main(int argc, char** argv)
{
    double x = 0;
    int command = 0;            // Введенная пользователем команда
    double **table = NULL;      // Таблицы

    printf("Выберите метод интерполяции: \n\
           1. сплайны \n\
           2. многомерная \n\
           > ");
    scanf("%d",&command);

    int code=0;
    switch (command)
    {
        case 1:
        {
            //input_limits(left,right,parts);

            table = (double**) malloc(sizeof(double*)*COLUMNS);
            
            for (int i=0; i<COLUMNS; i++)
                table[i] = (double*) malloc(sizeof(double)*SPLINE_SIZE);
            
            double x = 0;
            for (int i=0; i<SPLINE_SIZE; i++)
            {
                table[0][i] = x;
                table[1][i] = f(x);
                x += 1;
            }
            
            //fill_table(table,left,right,parts);
            
            puts("");
            print_table(table, SPLINE_SIZE, COLUMNS);
            puts("");
            
            printf("Введите x: ");
            scanf("%lf",&x);

            double res = spline(table[0], table[1], SPLINE_SIZE, x, code);
            if (code == 1)
            {
                printf("\n[Внимание] Выход за границы таблицы.\n");
                free(table);
                break;
            }
            
            printf("\nТочное значение:       %*.5f", 8, f(x));
            printf("\nВычисленное значение:  %*.5f\n\n", 8, res);
            
            free(table);
            break;
        }
        case 2:
        {
            double x = 0.0, y= 0.0;
            int nX = 0, nY = 0;
            
            // ****************************************
            // Заполение таблицы
            table = new (double* [SIZE+1]);
            
            for (int i=0; i<SIZE+1; i++)
                table[i] = new double[SIZE+1];
            
            table[0][0] = 0.0;
            
            for (int i=1; i<=SIZE; i++)
                table[0][i] = i;
            
            for (int i=1; i<=SIZE; i++)
                table[i][0] = i;
            
            for (int i=1; i<=SIZE; i++)
                for (int j=1; j<=SIZE; j++)
                    table[i][j] = table[i][0]*table[i][0] + table[0][j]*table[0][j];
            // ****************************************

            // ****************************************
            // Вывод таблицы
            puts("");
            puts("");
            for (int i=0; i<=SIZE; i++)
            {
                for (int j=0; j<=SIZE; j++)
                    if (i == 0 && j == 0)
                        printf("%8s |", "x/y");
                    else
                        printf("%8.2lf |",table[i][j]);
                
                
                printf("\n");
            }
            // ****************************************
            
            // ****************************************
            printf("\nВведите точку (x,y): ");
            scanf("%lf %lf", &x,&y);
            printf("\nВведите степень полинома для x: ");
            scanf("%d", &nX);
            printf("\nВведите степень полинома для y: ");
            scanf("%d", &nY);
            puts("");
            // ****************************************
            
            // ****************************************
            double res = multi(table, SIZE, SIZE, x, nX, y, nY, code);
            if (code == 1)
                printf("[Внимание] Произошла экстраполяция.\n");
            // ****************************************
            
            // ****************************************х
            printf("\nТочное значение для точки     (%.2lf, %.2lf): %*.3lf\n",   x, y, 8, x*x+y*y);
            printf("Вычислимое значение для точки   (%.2lf, %.2lf): %*.3lf\n\n", x, y, 8, res);
            // ****************************************

            free(table);
            break;
        }
    }


    return 0;
}
