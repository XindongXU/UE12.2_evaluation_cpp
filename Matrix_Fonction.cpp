#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <math.h>
// #include "Matrix.h"

class Matrix{
    public:
        std::tuple<int, int> shape();

        double Norme();
        double GetElement(int i, int j);
        void   SetElement(int i, int j, double x);
        void   show();
        void   initialise(int nb_row, int nb_col);

        double DivisScal(Matrix aux);
        Matrix MultiScal(double lambda);
        Matrix Somme(Matrix aux);
        Matrix Diffe(Matrix aux);
        Matrix Multi(Matrix aux);
        Matrix Transpose();


    protected:
        int nb_row;
        int nb_col;
        std::vector<double> elements;
};


void Matrix::initialise(int i, int j){ // pour initialiser une matrice de zeros de taille (i, j)
    nb_row = i;
    nb_col = j;
    elements.resize(nb_row * nb_col);
}

std::tuple<int, int> Matrix::shape(){ // cette méthode nous donne le shape de matrice;
    return {nb_row, nb_col};
}

double Matrix::GetElement(int i, int j){ // i = row dans math , j = col dans math mais pas celles dans l'informatique;
    return elements[nb_col * (i - 1) + (j - 1)];
}

void Matrix::SetElement(int i, int j, double x){ // i = row dans math , j = col dans math mais pas celles dans l'informatique;
    elements[nb_col * (i - 1) + (j - 1)] = x;
}

Matrix Matrix::Somme(Matrix aux){ // cette méthode réalise une opération de mat + aux;
    Matrix res;
    res.initialise(0, 0);
    //std::cout << (shape() == aux.shape()) << std::endl;
    if (shape() != aux.shape()){
        std::cout << "Shapes of matrix to operate are not propre, please check" << std::endl;
        return res;
    } else {
        res.initialise(nb_row, nb_col);
        // std::cout << nb_row << nb_col << std::endl;
        for(int i = 0; i < nb_row; i++){
            for(int j = 0; j < nb_col; j++){
                res.SetElement(i + 1, j + 1, GetElement(i + 1, j + 1) + aux.GetElement(i + 1, j + 1));
                // std::cout << i << j << std::endl;
            }
        }
        return res;
    }
}

Matrix Matrix::Diffe(Matrix aux){ // cette méthode réalise une opération de mat - aux;
    //std::cout << (shape() == aux.shape()) << std::endl;
    Matrix res;
    res.initialise(0, 0);
    if (shape() != aux.shape()){
        std::cout << "Shapes of matrix to operate are not propre, please check" << std::endl;
        return res;
    } else {
        res.initialise(nb_row, nb_col);
        // std::cout << nb_row << nb_col << std::endl;
        for(int i = 0; i < nb_row; i++){
            for(int j = 0; j < nb_col; j++){
                res.SetElement(i + 1, j + 1, GetElement(i + 1, j + 1) - aux.GetElement(i + 1, j + 1));
                // std::cout << i << j << std::endl;
            }
        }
        return res;
    }
}

Matrix Matrix::MultiScal(double lambda){ // cette méthode réalise une opération de lambda * mat;
    Matrix res;
    res.initialise(nb_row, nb_col);
    for (int i = 0; i < nb_row; i++){
        for (int j = 0; j < nb_col; j++){
            res.SetElement(i + 1, j + 1, GetElement(i + 1, j + 1)*lambda);
        }
    }
    return res;
}

Matrix Matrix::Multi(Matrix aux){ // cette méthode réalise une opération de mat * aux;
    auto m = nb_row; 
    auto n = nb_col;
    auto x = std::get<0>(aux.shape());
    auto y = std::get<1>(aux.shape());

    Matrix res;
    res.initialise(0, 0);
    if (n != x){
        std::cout << "Shapes of matrix to operate are not propre, please check" << std::endl;
        return res;
    } else {
        res.initialise(m, y);
        for (int i = 0; i < m; i++){
            for (int j = 0; j < y; j++){
                double e {0};
                for (int k = 0; k < n; k++){
                    e = e + GetElement(i + 1, k + 1) * aux.GetElement(k + 1, j + 1);
                }
                res.SetElement(i + 1, j + 1, e);
            }
        }
        return res;
    }
}

double Matrix::DivisScal(Matrix aux){ // cette méthode réalise une opération de division entre deux matrice qui ont tous la forme de (1, 1);
    double res = 0;
    if ((nb_col != 1)||(nb_row != 1)||(std::get<0>(aux.shape()) != 1)||(std::get<1>(aux.shape()) != 1)){
        std::cout << "Shapes of matrix to operate are not propre, please check" << std::endl;
        return res;
    }
    res = GetElement(1, 1)/aux.GetElement(1, 1);
    return res;
}

void Matrix::show(){ // cette méthode imprime la matrice;
    for(int i = 0 ; i < nb_row ; i++){
        for(int j = 0 ; j < nb_col ; j++){
            std::cout << GetElement(i + 1, j + 1) << " ";
            // std::cout << "yes" << " ";
        }
        std::cout << std::endl;
    }
}

double Matrix::Norme(){ // cette méthode réalise une opération de la norme euclidienne d'ordre 2;
    double res = 0;
    for (int i = 0; i < nb_row; i++){
        for (int j = 0; j < nb_col; j++){
            res = res + GetElement(i + 1, j + 1) * GetElement(i + 1, j + 1);
        }
    } 
    res = sqrt(res);
    return res;
}

Matrix Matrix::Transpose(){ // cette méthode réalise une opération de transpositon;
    Matrix res;
    res.initialise(nb_col, nb_row);
    for (int i = 0; i < nb_col; i++){
        for (int j = 0; j < nb_row; j++){
            res.SetElement(i + 1, j + 1, GetElement(j + 1, i + 1));
        }
    }
    return res;
}



Matrix condition(int N){ 
    // cette fonction nous donne la condition initiale de la distribution de chaleur, 
    // qui vérifie la relation de T(x, t) et les conditons aux limites;
    double dx {1.0/(N-1.0)};
    Matrix res;
    res.initialise(N, 1); // des matrices de {T}^(k) sont de forme de N lignes, 1 colonne;

    for (int i = 0; i < N; i++){
        res.SetElement(i + 1, 1, 1/2 + sin(2*M_PI*i*dx) - 1/2*cos(2*M_PI*i*dx));
        // relation de T(x ,t) dans l'énoncé;
    }
    res.SetElement(N, 1, 0.0); // conditions aux limites;
    return res;
}

void OutPut(int N, double dt, std::vector<Matrix> Tf){
    // cette fonction exporte des résultats calculés par la méthode de l'euler, dans un dossier data.txt
    // puis, en utilisant data_figures.ipynb, on trace des figures;
    std::ofstream dataout;
    dataout.open("data.txt");
    dataout << N << std::endl;
    dataout << dt << std::endl;
    for (int m = 0; m < Tf.size(); m++){
        for (int i = 0; i < std::get<0>(Tf[0].shape()); i++){
            auto ele = Tf[m].GetElement(i + 1, 1);
            dataout << ele << std::endl;
            }
        }
    dataout.close();
}

// int main(){ // on vérifie la classe et des méthodes avec des matrices données : mat et aux;
//     Matrix mat;
//     mat.initialise(2, 2);
//     mat.SetElement(1, 1, 6.0);
//     mat.SetElement(1, 2, 3.0);
//     mat.SetElement(2, 1, 2.0);
//     mat.SetElement(2, 2, 1.0);
//     // auto shape_mat = mat.shape();
//     std::cout << "mat :" << std::endl;
//     mat.show();
//     std::cout << "shape of mat = (" << std::get<0>(mat.shape()) << ", " << std::get<1>(mat.shape()) << ")" << std::endl;
//     std::cout << std::endl;

//     Matrix aux;
//     aux.initialise(2, 2);
//     aux.SetElement(1, 1, 4.0);
//     aux.SetElement(1, 2, 2.0);
//     aux.SetElement(2, 1, 1.0);
//     aux.SetElement(2, 2, 0.5);
//     std::cout << "aux :" << std::endl;
//     aux.show();
//     std::cout << std::endl;

//     std::cout << "operation of somme:" << std::endl;
//     Matrix res = mat.Somme(aux);
//     res.show();
//     std::cout << std::endl;

//     std::cout << "operation of difference:" << std::endl;
//     res = mat.Diffe(aux);
//     res.show();
//     std::cout << std::endl;

//     std::cout << "operation of multiplication with 0.5:" << std::endl;
//     res = mat.MultiScal(0.5);
//     res.show();
//     std::cout << std::endl;

//     std::cout << "operation of multiplication with aux:" << std::endl;
//     res = mat.Multi(aux);
//     res.show();
//     std::cout << std::endl;

//     std::cout << "operation of transposition:" << std::endl;
//     res = mat.Transpose();
//     res.show();
//     std::cout << std::endl;

//     std::cout << "operation of norm:" << std::endl;
//     std::cout << mat.Norme() << std::endl;


//     return 0;
// }