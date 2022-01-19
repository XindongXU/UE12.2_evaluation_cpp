#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include "Matrix_Fonction.cpp"
#include <math.h>

Matrix gradient_conjugue(Matrix A, Matrix b){
    // cette fonction nous aide résoudre des systèmes linéaires de forme de : A.X = b;
    Matrix X;
    X.initialise(std::get<0>(b.shape()), std::get<1>(b.shape()));
    Matrix r = b.Diffe(A.Multi(X));
    Matrix p = r;
    double beta;

    double alpha = r.Transpose().Multi(r).DivisScal(p.Transpose().Multi(A.Multi(p)));
    X = X.Somme(p.MultiScal(alpha));
    Matrix r1 = r.Diffe(A.Multi(p).MultiScal(alpha));
    
    while (r1.Norme() >= 0.001){
        beta = r1.Transpose().Multi(r1).DivisScal(r.Transpose().Multi(r));
        p = r1.Somme(p.MultiScal(beta));
        r = r1;
        alpha = r.Transpose().Multi(r).DivisScal(p.Transpose().Multi(A.Multi(p)));
        X = X.Somme(p.MultiScal(alpha));
        r1 = r.Diffe(A.Multi(p).MultiScal(alpha));
    }
    return X;
}

int main(){
    Matrix mat;
    mat.initialise(2, 2);
    mat.SetElement(1, 1, 2.0);
    mat.SetElement(1, 2, 3.0);
    mat.SetElement(2, 1, 6.0);
    mat.SetElement(2, 2, 7.0);
    std::cout << "A.x = B" << std::endl;
    std::cout << "A = " << std::endl; 
    mat.show();

    Matrix aux;
    aux.initialise(2, 1);
    aux.SetElement(1, 1, 5.0);
    aux.SetElement(2, 1, 13.);
    std::cout << "B = " << std::endl; 
    aux.show();

    Matrix res = gradient_conjugue(mat, aux);
    std::cout << "x = " << std::endl; 
    res.show();
    return EXIT_SUCCESS;
}