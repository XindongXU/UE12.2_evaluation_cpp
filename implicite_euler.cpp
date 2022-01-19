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
    int  N = 50;
    // on choisit N comme 50, pour garder la qualité de figure;
    auto Tt = condition(N);
    // on définit dans la matrice Tt la condition initiale {T}^(0);
    double dx {1.0/(N-1.0)};

    std::vector<Matrix> Tf;
    Tf.push_back(Tt);
    double tf = 0.5;
    double dt = 0.05;
    // quand on choisit dt, il faut vérifier que tf/dt soit même ordre de grandeur de N;
    // pour que la méthode de l'euler implicite fonctionne bien et que des résultats obtenus soient stables;
    Matrix K;
    K.initialise(N, N);
    // on alors définit la matrice K;
    // une remarque : il faut modifier la première et la dernière ligne pour que la température soit toujours nulle aux bords;
    // on suit exactement ce qui est écrit dans le sujet, 
    // sauf pour la première et dernière ligne : on mis -2 et que des zéros pour la première ligne, et des zéros puis un -2 sur la dernière;
    // de plus, il faut modifier la première et la dernière colonne pour que la matrice K soit toujours symétrique;
    for (int i = 0; i < N; i++){
        K.SetElement(i + 1, i + 1, -2./dx/dx);
    }
    for (int i = 1; i < N - 2; i++){
        K.SetElement(i + 2, i + 1, 1./dx/dx);
    }
    for (int i = 2; i < N - 1; i++){
        K.SetElement(i, i + 1, 1./dx/dx);
    }

    Matrix A;
    A.initialise(N, N);
    // on définit mtn une matrice A de 1 de la même forme de K;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            A.SetElement(i + 1, j + 1, 1.0);
        }
    }
    K = K.MultiScal(dt);
    A = A.Diffe(K);
    // A = (1 - dt * K)

    for (int i = 0; i < N-1; i++){
        A.SetElement(1, i+2, 0.0);
        A.SetElement(N, i+1, 0.0);
        A.SetElement(i+2, 1, 0.0);
        A.SetElement(i+1, N, 0.0);
    }
    // A.show();
    
    // Euler implicite
    for(double t = 0; t < tf; t = t + dt){
        Tt = gradient_conjugue(A, Tt);
        // (1 - dt * K).{T}^(k+1) = {T}^(k)
        // A.{T}^(k+1) = {T}^(k)
        // Tt.SetElement(1, 1, 0.0);
        // Tt.SetElement(N, 1, 0.0);
        // on peut rassurer que les conditions aux limites sont atteintes, mais pas trop nécessaire grâce à notre bonne définition de K;
        Tf.push_back(Tt);
    }
    int zhubao = (tf/dt);
    // on affiche la distribution de température finale;
    std::cout << zhubao << std::endl;
    Tf[zhubao].show();
    
    OutPut(N, dt, Tf);

    return EXIT_SUCCESS;
}