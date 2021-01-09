#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

class Interval{
    public:
        double from;
        double to;

        Interval(double from, double to){
            this->from = from;
            this->to = to;
        }
};

class EiFunction{     // y(x) = a*x + b
    private:
        int i;
        int n;

    public:
        EiFunction(int i, int n){
            this->i = i;
            this->n = n;
        }

        double normal(double x){
            return (x >= 2.0*(i-1)/n && x <= 2.0*i/n) ? n*0.5*x-i+1 : ( (x >= 2.0*i/n && x <= 2.0*(i+1)/n) ? -n*0.5*x+i+1 : 0);
        }

        double derivative(double x){      // pochodna funkcji
            return (x >= 2.0*(i-1)/n && x <= 2.0*i/n) ? n*0.5 : ( (x >= 2.0*i/n && x <= 2.0*(i+1)/n) ? -n*0.5 : 0);
        }
};

const Interval omega(0,2);      // Przedział w którym szukamy funkcji
const int integralAcc = 50000;       // Dokładność obliczania całki

int k(double x){
    if(x >= 0 && x <= 1)    return 1;
    if(x > 1 && x <= 2) return 2;
    return 0;
}

double calculateIntegral(EiFunction u, EiFunction v){       // Metoda kwadratów
    double score = 0;
    double h = (omega.to - omega.from)/(double) integralAcc;
    double x;

    for(int i = 1; i <= integralAcc; i++){
        x = omega.from + i*h;
        // score += k(x)*u.derivative(x)*v.derivative(x)*h;
        score += u.derivative(x)*v.derivative(x)*h;

    }

    return score;
}

double B(EiFunction u, EiFunction v){
    return calculateIntegral(u, v) - (u.normal(0)*v.normal(0));
}

double L(EiFunction v){
    return -20*v.normal(0);
}


double **initialize2DMatrix(int n){
    double **A = new double*[n];
    
    for(int i = 0; i < n; i++)
        A[i] = new double [n+1];

    return A;
}

void delete2DMatrix(double **A, int n){
    for(int i = 0; i < n; i++)
        delete [] A[i];

    delete [] A;
}

void print2DMatrix(double **A, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n+1; j++){
            cout.width( 10 );
            cout.fill( ' ' );
            cout << A[i][j];
        }
        cout << endl;
    }
}

void print1DMatrix(double *tab, int n){
    for(int i = 0; i < n; i++){
            cout.width( 10 );
            cout.fill( ' ' );
            cout << tab[i];
    }
    cout << endl;
}

double ** makingMatrix(int n){
    double **A = initialize2DMatrix(n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)
            A[i][j] = calculateIntegral(EiFunction(i, n), EiFunction(j, n));
        
        A[i][n] = L(EiFunction(i, n));
    }

    return A;
}

double *solve2DMatrix(double **A, int n){       // Na podstawie algorytmu eliminacji Gaussa
    double *sollution = new double[n];
    double tmp;

    for(int j = 0; j < n-1; j++){
        for(int i = j+1; i < n; i++){
            tmp = A[i][j]/A[j][j];

            for(int k =0; k < n+1; k++){
                A[i][k] -= A[j][k]*tmp;
            }
        }
    }

    for(int i = n-1; i >= 0; i--){
        double s = 0;
        for(int j = i+1; j < n; j++)
            s += A[i][j] * sollution[j];
        
        sollution[i] = (A[i][n] - s)/A[i][i];
    }
    
    return sollution;
}

double calculateXFromSollutionTab(double *sollution, int n, double x){
    double score = 0;

    for(int i = 0; i < n; i++){
        score += sollution[i]*EiFunction(i,n).normal(x);
    }

    return score;
}

void makeFileToPlot(double *sollution, int n){
    ofstream file("data.txt");

    double h = (omega.to - omega.from)/(double) integralAcc;
    double x;

    for(int i = 1; i <= integralAcc; i++){
        x = omega.from + i*h;
        file.width(10);
        file.fill(' ');
        file << x;
        file.width(10);
        file.fill(' ');
        file << calculateXFromSollutionTab(sollution, n, x) << endl;
    }

    file.close();
}

void calculatingSolution(int n){
    double **A = makingMatrix(n);
    print2DMatrix(A, n);
    cout << endl;
    double *sollution = solve2DMatrix(A, n);

    print2DMatrix(A, n);
    cout << endl;
    print1DMatrix(sollution, n);

    makeFileToPlot(sollution, n);

    delete [] sollution;
    delete2DMatrix(A,n);
}


int main(){
    int n;
    cout << "Podaj ilość funkcji bazowych: ";
    cin >> n;
    // EiFunction e1(1, n);
    // EiFunction e2(2, n);
    // cout << B(e1, e1) << endl;
    // cout << B(e1, e2) << endl;
    // cout << B(e2, e1) << endl;
    // cout << B(e2, e2) << endl;
    cout << (double)0*(-1) << endl;

    calculatingSolution(n);

    return 0;
}