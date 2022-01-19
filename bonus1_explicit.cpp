#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include "Matrix_Fonction.cpp"
#include <math.h>
#include <string>
#include <time.h>

template<typename T>
T RandomGenerate(T nbmin, T nbmax)
{
	T temp;
	if (nbmin > nbmax)
	{
		temp = nbmax;
		nbmax = nbmin;
		nbmin = temp;
	}
	return rand() / (double)RAND_MAX *(nbmax - nbmin) + nbmin;
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
    double dt = 0.5/N/100;
    // quand on choisit dt, il faut vérifier que tf/dt soit beaucoup plus grand que N;
    // pour que la méthode de l'euler explicite fonctionne bien et que des résultats obtenus soient stables;
    
    std::vector<double> D(N);
    for (int i = 0; i < N; i++){
        D[i] = RandomGenerate<double>(0.5, 1.5);
    }

    Matrix K;
    K.initialise(N, N);
    for (int i = 0; i < N; i++){
        K.SetElement(i + 1, i + 1, (-D[i+1]-D[i+2])/dx/dx);
    }
    for (int i = 1; i < N - 2; i++){
        K.SetElement(i + 2, i + 1, D[i+2]/dx/dx);
    }
    for (int i = 2; i < N - 1; i++){
        K.SetElement(i, i + 1, D[i+1]/dx/dx);
    }

    // Euler explicite
    for(double t = 0; t < tf; t = t + dt){
        Tt = Tt.Somme((K.Multi(Tt)).MultiScal(dt));
        Tf.push_back(Tt);
    }

    int zhubao = (tf/dt-1);
    Tf[zhubao].show();
    // on affiche la distribution de température finale;

    OutPut(N, dt, Tf);
    return EXIT_SUCCESS;
}