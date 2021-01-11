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

const Interval omega(0,2);      // Przedział w którym szukamy rozwiązania
const int pointsAcc = 50000;       // Ilość punktów do wyrysowania wykresu

const double eps = 1e-12; // stała przybliżenia zera

class EiFunction{     // y(x) = a*x + b 
    public:
        int i;      // numer funkcji
        int n;

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

int k(double x){
    if(x >= 0 && x <= 1)    return 1;
    if(x > 1 && x <= 2) return 2;
    return 0;
}

void setABforIntegral(double &a, double &b, EiFunction u, EiFunction v){        // ustawia przedziały całkowania
    if(u.i == v.i){
        a = max(omega.from, 2.0*(u.i-1)/u.n);
        b = min(omega.to, 2.0*(u.i+1)/u.n);
    } else if (u.i + 1 == v.i){     // u.i <= v.i
        a = 2.0*(v.i-1)/v.n;
        b = 2.0*v.i/v.n;
    } else {
        a = 0;
        b = 0;
    }
}

// Oblicza wartość całki dla funkcji u i v dla przedziału [a, b]
double calculateIntegralFromAtoB(EiFunction u, EiFunction v, double a, double b){
    double x = 5.77350269189625764507e-01;
    double w = 1.0;

    double c = 0.5 * (b-a);
    double d = 0.5 * (b+a);
    double dum = c * x;

    return c * w * ( (k(d-dum)*u.derivative(d-dum)*v.derivative(d-dum)) + (k(d+dum)*u.derivative(d+dum)*v.derivative(d+dum)) );
}

// Oblicza wartość całki dla funkcji u, v
double calculateIntegral(EiFunction u, EiFunction v){
    double a, b;
    setABforIntegral(a,b,u,v);

    if(a == 0 && b == 0)    return 0;
    else if(a <= 1.0 && b > 1.0 )    return calculateIntegralFromAtoB(u,v,a,1.0) + calculateIntegralFromAtoB(u,v,1.0,b);
    else return calculateIntegralFromAtoB(u,v,a,b);
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
            cout.width( 12 );
            cout.fill( ' ' );
            cout << A[i][j];
        }
        cout << endl;
    }
}

void print1DMatrix(double *tab, int n){
    for(int i = 0; i < n; i++){
            cout.width( 8 );
            cout.fill( ' ' );
            cout << tab[i];
    }
    cout << endl;
}

// Tworzy i wypełnia macierz
double ** makingMatrix(int n){
    double **A = initialize2DMatrix(n);

    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++)
            A[j][i] = A[i][j] = B(EiFunction(i, n), EiFunction(j, n));
        
        A[i][n] = L(EiFunction(i, n));
    }

    return A;
}

// Oblicza rozwiązania na podstawie algorytmu eliminacji Gaussa
double *solve2DMatrix(double **A, int n){
    double *sollution = new double[n];

    for(int i = 0; i < n; i++)
        sollution[i] = 0;

    double tmp;

    for(int i = 0; i < n; i++){
        for(int k = i+1; k < n; k++){
            if(fabs(A[i][i]) <= fabs(A[k][i])){
                for(int j = 0; j<=n; j++){
                    swap(A[i][j], A[k][j]);
                }
            }
        }
    }

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

// Oblicza wartość x na podstawie wyliczonego rozwiązania
double calculateXFromSollutionTab(double *sollution, int n, double x){
    double score = 0;

    for(int i = 0; i < n; i++){
        score += sollution[i]*EiFunction(i,n).normal(x);
    }

    return score;
}

// Tworzy plik z punktami do wyrysowania przez gnuplot
void makeFileToPlot(double *sollution, int n){
    ofstream file("data.txt");
    file << "plot '-'\n";

    double h = (omega.to - omega.from)/(double) pointsAcc;
    double x;

    for(int i = 1; i <= pointsAcc; i++){
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

    calculatingSolution(n);
    system("gnuplot -p < \"data.txt\"");

    return 0;
}