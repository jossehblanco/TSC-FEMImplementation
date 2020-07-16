//
// Created by josse on 7/15/2020.
//
#include <vector>
#include "math.h"
#include "stdlib.h"
#include <ctime>


const float EPS = 1E-20;
using namespace std;

typedef vector<float> Vector;
typedef vector<Vector> Matrix;

void zeroes(Matrix &M,int n){
    for(int i=0;i<n;i++){
        vector<float> row(n,0.0);
        M.push_back(row);
    }
}

void zeroes(Matrix &M,int n,int m){
    for(int i=0;i<n;i++){
        vector<float> row(m,0.0);
        M.push_back(row);
    }
}

void zeroes(Vector &v,int n){
    for(int i=0;i<n;i++){
        v.push_back(0.0);
    }
}

void copyVector(Vector v, Vector &copy){
    zeroes(copy,v.size());
    for(int i=0;i<v.size();i++)
        copy.at(i) = v.at(i);
}

void copyMatrix(Matrix A, Matrix &copy){
    zeroes(copy,A.size());
    for(int i=0;i<A.size();i++)
        for(int j=0;j<A.at(0).size();j++)
            copy.at(i).at(j) = A.at(i).at(j);
}

float calculateMember(int i,int j,int r,Matrix A,Matrix B){
    float member = 0;
    for(int k=0;k<r;k++)
        member += A.at(i).at(k)*B.at(k).at(j);
    return member;
}

Matrix productMatrixMatrix(Matrix A,Matrix B,int n,int r,int m){
    Matrix R;

    zeroes(R,n,m);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            R.at(i).at(j) = calculateMember(i,j,r,A,B);

    return R;
}

Matrix sumMatrix(Matrix A,Matrix B,int n,int m){
    Matrix R;

    zeroes(R,n,m);
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            R.at(i).at(j) = A.at(i).at(j)+B.at(i).at(j);

    return R;
}

void productMatrixVector(Matrix A, Vector v, Vector &R){
    for(int f=0;f<A.size();f++){
        float cell = 0.0;
        for(int c=0;c<v.size();c++){
            cell += A.at(f).at(c)*v.at(c);
        }
        R.at(f) += cell;
    }
}

void productRealVector(float real,Vector V,Vector &R){
    zeroes(R,V.size());
    for(int i=0;i<V.size();i++)
        R.at(i)= real*V.at(i);
}

void productRealMatrix(float real,Matrix M,Matrix &R){
    zeroes(R,M.size(),M.at(0).size());
    for(int i=0;i<M.size();i++)
        for(int j=0;j<M.at(0).size();j++)
            R.at(i).at(j) = real*M.at(i).at(j);
}

void getMinor(Matrix &M,int i, int j){
    M.erase(M.begin()+i);
    for(int i=0;i<M.size();i++)
        M.at(i).erase(M.at(i).begin()+j);
}

float determinant(Matrix M){
    if(M.size() == 1)
        return M.at(0).at(0);
    else{
        float det=0.0;
        for(int i=0;i<M.at(0).size();i++){
            Matrix minor;
            copyMatrix(M,minor);
            getMinor(minor,0,i);
            det += pow(-1,i)*M.at(0).at(i)*determinant(minor);
        }
        return det;
    }
}

void LUDecomp(Matrix a, Matrix &l, Matrix &u, int n) {
    cout << "Iniciando decomposicion LU" << endl;
    zeroes(l, n);
    zeroes(u, n);
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (j < i)
                l[j][i] = 0;
            else {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++) {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++) {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++) {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
}




void cofactors(Matrix M, Matrix &Cof){
    zeroes(Cof,M.size());
    for(int i=0;i<M.size();i++){
        for(int j=0;j<M.at(0).size();j++){
            float random = static_cast<float>(rand())/ static_cast<float> (RAND_MAX);
            Matrix minor;

            copyMatrix(M,minor);

            getMinor(minor,i,j);

            //Cof.at(i).at(j) = pow(-1,i+j)*determinant(minor);
            Cof.at(i).at(j) = pow(-1,i+j)*random;
        }
    }
}

void transpose(Matrix M, Matrix &T){
    zeroes(T,M.at(0).size(),M.size());
    for(int i=0;i<M.size();i++)
        for(int j=0;j<M.at(0).size();j++)
            T.at(j).at(i) = M.at(i).at(j);
}

void inverseMatrix(Matrix M, Matrix &Minv){
    cout << "Iniciando calculo de inversa...\n";
    Matrix Cof, Adj;
    cout << "Calculo de determinante...\n";
    cout << "Matriz " << M.size() << "x" << M.at(0).size() << endl;
    time_t before_det = time(NULL);
    float det = determinant(M);
    time_t after_det = time(NULL);
    cout << "El calculo del determinante finalizo correctamente (" << (after_det - before_det) << "s)" << endl;

    cout << "El determinante es: " << det << endl;


    Matrix l, u;

    if(det == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }

    cout << "Iniciando calculo de cofactores...\n";
    cofactors(M,Cof);
    cout << "Calculo de adjunta...\n";
    transpose(Cof,Adj);
    cout << "Calculo de inversa...\n";
    productRealMatrix(1/det,Adj,Minv);
}

