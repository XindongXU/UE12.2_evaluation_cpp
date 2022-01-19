#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
// #include "Matrix.h"
#include "Matrix_Fonction.cpp"
#include <math.h>

int main(){
    int  N = 50;
    // on choisit N comme 50, pour garder la qualité de figure;
    auto Tt = condition(N);
    // on définit dans la matrice Tt la condition initiale {T}^(0);
    double dx {1.0/(N-1.0)};

    std::vector<Matrix> Tf;
    Tf.push_back(Tt);
    double tf = 0.5;
    double dt = 0.5/N/100;
    // quand on choisit dt, il faut vérifier que tf/dt soit beaucoup plus grand que N;
    // pour que la méthode de l'euler explicite fonctionne bien et que des résultats obtenus soient stables;
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

    // Euler explicite
    for(double t = 0; t < tf; t = t + dt){
        Tt = Tt.Somme((K.Multi(Tt)).MultiScal(dt));
        // {T}^(k+1) = {T}^(k) - dt * K.{T}^(k);

        // Tt.SetElement(1, 1, 0.0);
        // Tt.SetElement(N, 1, 0.0);
        // on peut rassurer que les conditions aux limites sont atteintes, mais pas trop nécessaire grâce à notre bonne définition de K;
        Tf.push_back(Tt);
    }

    int zhubao = (tf/dt/2);
    std::cout << zhubao << std::endl;
    Tf[zhubao].show();
    // on affiche la distribution de température finale;

    OutPut(N, dt, Tf);
    return EXIT_SUCCESS;
}