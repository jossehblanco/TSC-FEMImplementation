//
// Created by josse on 7/15/2020.
//

#define _USE_MATH_DEFINES
#include <math.h>

float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);

    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}


void calculateAlphaPrime(int i,Matrix &Alpha_prime,mesh m){
    zeroes(Alpha_prime,3);
    element e = m.getElement(i);

    Alpha_prime.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    Alpha_prime.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    Alpha_prime.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    Alpha_prime.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    Alpha_prime.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    Alpha_prime.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    Alpha_prime.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    Alpha_prime.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    Alpha_prime.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateOmegaPrime(Matrix &OmegaP){
    zeroes(OmegaP,3,12);

    OmegaP.at(0).at(0) = -1;
    OmegaP.at(0).at(1) =  1;
    OmegaP.at(0).at(2) =  0;
    OmegaP.at(0).at(3) =  0;
    OmegaP.at(0).at(4) = -1;
    OmegaP.at(0).at(5) =  1;
    OmegaP.at(0).at(6) =  0;
    OmegaP.at(0).at(7) =  0;
    OmegaP.at(0).at(8) = -1;
    OmegaP.at(0).at(9) =  1;
    OmegaP.at(0).at(10) =  0;
    OmegaP.at(0).at(11) =  0;

    OmegaP.at(1).at(0) = -1;
    OmegaP.at(1).at(1) = 0;
    OmegaP.at(1).at(2) = 1;
    OmegaP.at(1).at(3) = 0;
    OmegaP.at(1).at(4) = -1;
    OmegaP.at(1).at(5) = 0;
    OmegaP.at(1).at(6) = 1;
    OmegaP.at(1).at(7) = 0;
    OmegaP.at(1).at(8) = -1;
    OmegaP.at(1).at(9) = 0;
    OmegaP.at(1).at(10) = 1;
    OmegaP.at(1).at(11) = 0;

    OmegaP.at(2).at(0) = -1;
    OmegaP.at(2).at(1) = 0;
    OmegaP.at(2).at(2) = 0;
    OmegaP.at(2).at(3) = 1;
    OmegaP.at(2).at(4) = -1;
    OmegaP.at(2).at(5) = 0;
    OmegaP.at(2).at(6) = 0;
    OmegaP.at(2).at(7) = 1;
    OmegaP.at(2).at(8) = -1;
    OmegaP.at(2).at(9) = 0;
    OmegaP.at(2).at(10) = 0;
    OmegaP.at(2).at(11) = 1;
}

void calculateBeta(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1;
    C.at(0).at(1) = 1;
    C.at(0).at(2) = 0;
    C.at(0).at(3) = 0;

    C.at(1).at(0) = -1;
    C.at(1).at(1) = 0;
    C.at(1).at(2) = 1;
    C.at(1).at(3) = 0;

    C.at(2).at(0) = -1;
    C.at(2).at(1) = 0;
    C.at(2).at(2) = 0;
    C.at(2).at(3) = 1;
}


void calculateGamma(Matrix &m){
    zeroes(m,12,3);

    m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
    m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0;
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
    m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0;
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
    m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0;
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
    m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0;
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
    m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1;
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
    m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1;
}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);

    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateSigma(Vector &s){
    zeroes(s, 4);
    s.at(0) = 1;
    s.at(1) = 1;
    s.at(2) = 1;
    s.at(3) = 1;

}


void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(M_X);
    f.at(1) += m.getParameter(M_Y);
    f.at(2) += m.getParameter(M_Z);

}


Matrix createLocalM(int e,mesh &m){
    Matrix matrixK,matrixU,matrixC,matrixB, matrixD, matrixE;

    float u_bar,nu,rho,Ve,J,Determinant;

    /*Se tiene el siguiente sistema perteneciente a la matriz M
     * [ -K+U  -C+B ]
       [   D     E  ]
    */


    //Definiendo terminos comunes
    Matrix  gamma_ones, gamma_full, alphaprime, omegaprime, alphaprime_trans, omegap_trans, beta, beta_trans;
    Vector sigma, f, f_final;
    calculateGamma(gamma_ones);
    float gamma_real = calculateLocalJ(e, m)/24;
    productRealMatrix(gamma_real, gamma_ones, gamma_full);
    calculateAlphaPrime(e, alphaprime, m);
    calculateOmegaPrime(omegaprime);
    calculateSigma(sigma);
    transpose(alphaprime, alphaprime_trans);
    transpose(omegaprime, omegap_trans);
    calculateBeta(beta);
    transpose(beta, beta_trans);
    calculateF(f, m);
    float determinante = calculateLocalD(e, m);
    float volumen = calculateLocalVolume(e, m);
    float jacobiano = calculateLocalJ(e, m);


    //Saliendo si el J = 0

    if(jacobiano == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }



    //Definiendo matriz K
    //12x12
    // k*V/d^2 omega_trans alphaprime_trans, alphaprime, omega
    float matrixK_real = (-1) * (m.getParameter(KAPPA) * volumen)/ (determinante * determinante);
    productRealMatrix(matrixK_real,productMatrixMatrix(omegap_trans,  productMatrixMatrix(alphaprime_trans, productMatrixMatrix(alphaprime,  omegaprime, 3,3,  12), 3, 3, 12 ), 12,3,  12), matrixK);

    //Definiendo matriz U
    //(2*alpha*J)/d

    float matrixU_real = (2*m.getParameter(ALPHA)*jacobiano)/determinante;
    productRealMatrix(matrixU_real, productMatrixMatrix(gamma_full, productMatrixMatrix(alphaprime, omegaprime, 3, 3, 12), 12, 3, 12), matrixU);
    //12x12


    //Definiendo matriz C
    float matrixC_real = (-1) * ( m.getParameter(CCONST)*m.getParameter(CCONST)* M_PI * volumen )/ (pow(M_E, m.getParameter(OMEGA))* determinante * determinante);
    productRealMatrix(matrixC_real, productMatrixMatrix(omegap_trans, productMatrixMatrix(alphaprime_trans, productMatrixMatrix(alphaprime, beta, 3, 3 ,4), 3, 3, 4), 12, 3, 4), matrixC);
    //12x4

    //Definiendo matriz B
    float matrixB_real = (m.getParameter(OMEGA) * jacobiano)/( pow(cos(m.getParameter(PHI)), 5)* determinante);
    productRealMatrix(matrixB_real, productMatrixMatrix(gamma_full, productMatrixMatrix(alphaprime, beta, 3, 3, 4), 12, 3, 4), matrixB);
    //12x4


    //Definiendo Matriz D

    float matrixD_real = (m.getParameter(ALPHA) * volumen)/(determinante * determinante);
    productRealMatrix(matrixD_real, productMatrixMatrix(beta_trans, productMatrixMatrix(alphaprime_trans, productMatrixMatrix(alphaprime, omegaprime, 3, 3, 12), 3, 3, 12), 4, 3, 12), matrixD);
    //4x12

    //Definiendo Matriz E
    float matrixE_real = (-1)* (log(M_PI) * volumen)/(determinante * determinante);
    productRealMatrix(matrixE_real, productMatrixMatrix(beta_trans, productMatrixMatrix(alphaprime_trans, productMatrixMatrix(alphaprime, beta, 3, 3, 4), 3, 3, 4), 4, 3, 4), matrixE);
    //4x4

    //Matriz M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixK,matrixU,12,12));
    ubicarSubMatriz(M,0,11,12,15, sumMatrix(matrixC, matrixB, 12, 4));
    ubicarSubMatriz(M,12,15,0,11,matrixD);
    ubicarSubMatriz(M,12,15,12,15,matrixE);
    return M;
}



Vector createLocalb(int e,mesh &m){
    float jacobiano;
    Vector b_final,b,f_vect;
    Matrix gamma_ones;
    calculateGamma(gamma_ones);
    calculateF(f_vect, m);

    jacobiano = calculateLocalJ(e,m);

    if(jacobiano == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }

    zeroes(b,16);
    productMatrixVector(gamma_ones,f_vect,b);
    //16x1
    productRealVector(jacobiano/24,b,b_final);

    return b_final;
}

