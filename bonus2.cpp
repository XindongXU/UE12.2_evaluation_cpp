#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include "Matrix_Fonction.cpp"
#include <math.h>
#include <chrono>   
using namespace std;
using namespace chrono;

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
    // Euler implicite
    int  N = 50;
    auto Tt = condition(N);
    double dx {1.0/(N-1.0)};
    std::vector<Matrix> Tf;
    Tf.push_back(Tt);
    double tf = 0.5;
    double dt = 0.05;
    Matrix K;
    K.initialise(N, N);
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

    auto start = system_clock::now();
    for(double t = 0; t < tf; t = t + dt){
        Tt = gradient_conjugue(A, Tt);
        Tf.push_back(Tt);
    }

    auto end   = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout <<  "Méthode de l'euler implicite, ça prend " << double(duration.count()) * microseconds::period::num / microseconds::period::den << "seconds" << endl;
    
////////////////////////////////////////////////

    // Euler implicite
    N = 50;
    Tt = condition(N);
    dx = 1.0/(N-1.0);

    std::vector<Matrix> Tf2;
    Tf2.push_back(Tt);
    tf = 0.5;
    dt = 0.5/N/100;
    for (int i = 0; i < N; i++){
        K.SetElement(i + 1, i + 1, -2./dx/dx);
    }
    for (int i = 1; i < N - 2; i++){
        K.SetElement(i + 2, i + 1, 1./dx/dx);
    }
    for (int i = 2; i < N - 1; i++){
        K.SetElement(i, i + 1, 1./dx/dx);
    }

    start = system_clock::now();
    for(double t = 0; t < tf; t = t + dt){
        Tt = Tt.Somme((K.Multi(Tt)).MultiScal(dt));
        Tf2.push_back(Tt);
    }

    end   = system_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout <<  "Méthode de l'euler explicite, ça prend " << double(duration.count()) * microseconds::period::num / microseconds::period::den << "seconds" << endl;

    return EXIT_SUCCESS;
}