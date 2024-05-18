#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    // ToDo
    if(a.rows!=b.rows || a.cols!=b.cols){
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }else{
        Matrix c;
        c.rows = a.rows;
        c.cols = a.cols;
        for(int i = 0; i < c.rows; i++){
            for(int j =0; j < c.cols; j++){
                c.data[i][j] = a.data[i][j] + b.data[i][j]; 
            }
        }
        return c;
    }
    
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    // ToDo
    if(a.rows!=b.rows || a.cols!=b.cols){
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }else{
        Matrix c;
        c.rows = a.rows;
        c.cols = a.cols;
        for(int i = 0; i < c.rows; i++){
            for(int j =0; j < c.cols; j++){
                c.data[i][j] = a.data[i][j] - b.data[i][j]; 
            }
        }
        return c;
    }
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    // ToDo
    if(a.cols != b.rows){
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }else{
        Matrix c;
        c.rows = a.rows;
        c.cols = b.cols;
        for(int i = 0; i < c.rows; i++){
            for(int j = 0; j < c.cols; j++){
                c.data[i][j] = 0.0;
                for(int k = 0; k < c.cols; k++){
                    c.data[i][j] += a.data[i][k] + b.data[i][k];
                }
            }
        }
        return c;
    }
    
}

Matrix scale_matrix(Matrix a, double k)
{
    // ToDo
    Matrix c;
    c.rows = a.rows;
    c.cols = a.cols;
    for(int i = 0; i < a.rows; i++){
        for(int j = 0; j < a.cols; j++){
            c.data[i][j] = k * a.data[i][j];
        }
    }
    return c;
}

Matrix transpose_matrix(Matrix a)
{
    // ToDo
    Matrix c;
    c.rows = a.cols;
    c.cols = a.rows;
    for(int i = 0; i < c.rows; i++){
        for(int j = 0; j < c.cols; j++){
            c.data[i][j] = a.data[j][i];
        }
    }
    return c;
}

double det_matrix(Matrix a)
{
    // ToDo
    if(a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }else{
        if(a.rows == 2){
            return (a.data[0][0]*a.data[1][1] - a.data[0][1]*a.data[1][0]);
        }else{
            double det;
            for(int j = 0; j < a.cols; j++){
                Matrix a0j;
                a0j.rows = a.rows - 1;
                a0j.cols = a.cols - 1;
                for(int i0 = 1; i0 < a.rows; i0++){
                    for(int j0 = 0; j0 < a.cols; j0++){
                        if(j0 < j){
                            a0j.data[i0-1][j0] = a.data[i0][j0];
                        }else if(j0 == j){
                            continue;
                        }else{
                            a0j.data[i0-1][j0-1] = a.data[i0][j0];
                        }
                    }
                }
                det += a.data[0][j] * (j % 2 ? 1 : -1) * det_matrix(a0j);
            }
            return det;
        }
    }
    return 0;
}

Matrix inv_matrix(Matrix a)
{
    // ToDo
    if(a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }else if(det_matrix(a) == 0){
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }else{
        Matrix c;
        c.rows = a.rows;
        c.cols = a.cols;
        for(int i = 0; i < c.rows; i++){
            for(int j = 0; j < c.cols; j++){
                Matrix aij;
                aij.rows = a.rows - 1;
                aij.cols = a.cols - 1;
                for(int i0 = 0; i0 < a.rows; i0++){
                    for(int j0 = 0; j0 < a.cols; j0++){
                        if(i0 < i && j0 < j){
                            aij.data[i0][j0] = a.data[i0][j0];
                        }else if(i0 < i && j0 > j){
                            aij.data[i0][j0-1] = a.data[i0][j0];
                        }else if(i0 > i && j0 < j){
                            aij.data[i0-1][j0] = a.data[i0][j0];
                        }else if(i0 > i && j0 > j){
                            aij.data[i0-1][j0-1] = a.data[i0][j0];
                        }else{
                            continue;
                        }
                    }
                }
                c.data[i][j] = ((i + j) % 2 ? 1 : -1) * det_matrix(aij)/det_matrix(a);
            }
        }
        return c;
    }
    
}

int rank_matrix(Matrix a)
{
    // ToDo
    int r;
    r = ((a.rows >= a.cols) ? a.cols : a.rows);
    if(a.rows > a.cols){
        a = transpose_matrix(a);
    }
    for(int i = 0; i < a.rows; i++){
        if(a.data[i][i] == 0){  //若aii==0，向下寻找第i列不为0的行
            int k = i+1;
            while(a.data[k][i] == 0 && k < a.rows){
                k++;
            }
            for(int j = 0; j < a.cols; j++){  //找到后交换第i行与第k行
                double temp;
                temp = a.data[i][j];
                a.data[i][j] = a.data[k][j];
                a.data[k][j] = temp;
            }
        }
    }
    for(int i = 0; i < a.rows; i++){   //从第0行开始进行高斯消元
        for(int k = i + 1; k < a.rows; k++){  // 对于当前主对角线上的元素aii，将其下方的所有元素通过行运算消除为零，使得当前列下方的元素全部为零。
            for(int j = 0; j < a.cols; j++){
                a.data[k][j] = a.data[k][j] - (a.data[k][i]/a.data[i][i]) * a.data[i][j];
            }
        }
        if(a.data[i + 1][i + 1] == 0){   //若下一行主对角线上元素为0，向下寻找这一列元素不为0的行，并与下一行交换，以便进行下一轮高斯消元，这是下一轮消元前的准备工作
            int k = i + 2;
            while(a.data[k][i] == 0 && k < a.rows){
                k++;
            }
            if(k == a.rows){    //若这一列下方的元素全部为0，需把这一列交换至最后一列，并将秩减一
                for(int i0 = 0; i0 < a.rows; i0++){
                    double temp = a.data[i0][i + 1];
                    for(int j0 = i + 1; j0 < a.cols - 1; j0++){
                        a.data[i0][j0] = a.data[i0][j0 + 1];
                    }
                    a.data[i0][a.cols-1] = temp; 
                }
                r--;
            }else{
                for(int j = 0; j < a.cols; j++){  //若能找到非零元素，与下一行交换
                    double temp;
                    temp = a.data[i + 1][j];
                    a.data[i + 1][j] = a.data[k][j];
                    a.data[k][j] = temp;
                }
            }
        }
    }
    return r;
}

double trace_matrix(Matrix a)
{
    // ToDo
    if(a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }else{
        double tr = 0.0;
        for(int i = 0; i < a.cols; i++){
            tr += a.data[i][i];
        }
        return tr;
    }
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}