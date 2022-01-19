#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include "Matrix_Fonction.cpp"
#include <math.h>

int main(){ // on vérifie la classe et des méthodes avec des matrices données : mat et aux;
    Matrix mat;
    mat.initialise(2, 2);
    mat.SetElement(1, 1, 6.0);
    mat.SetElement(1, 2, 3.0);
    mat.SetElement(2, 1, 2.0);
    mat.SetElement(2, 2, 1.0);
    // auto shape_mat = mat.shape();
    std::cout << "mat :" << std::endl;
    mat.show();
    std::cout << "shape of mat = (" << std::get<0>(mat.shape()) << ", " << std::get<1>(mat.shape()) << ")" << std::endl;
    std::cout << std::endl;

    Matrix aux;
    aux.initialise(2, 2);
    aux.SetElement(1, 1, 4.0);
    aux.SetElement(1, 2, 2.0);
    aux.SetElement(2, 1, 1.0);
    aux.SetElement(2, 2, 0.5);
    std::cout << "aux :" << std::endl;
    aux.show();
    std::cout << std::endl;

    std::cout << "operation of somme:" << std::endl;
    Matrix res = mat.Somme(aux);
    res.show();
    std::cout << std::endl;

    std::cout << "operation of difference:" << std::endl;
    res = mat.Diffe(aux);
    res.show();
    std::cout << std::endl;

    std::cout << "operation of multiplication with 0.5:" << std::endl;
    res = mat.MultiScal(0.5);
    res.show();
    std::cout << std::endl;

    std::cout << "operation of multiplication with aux:" << std::endl;
    res = mat.Multi(aux);
    res.show();
    std::cout << std::endl;

    std::cout << "operation of transposition:" << std::endl;
    res = mat.Transpose();
    res.show();
    std::cout << std::endl;

    std::cout << "operation of norm:" << std::endl;
    std::cout << mat.Norme() << std::endl;


    return 0;
}