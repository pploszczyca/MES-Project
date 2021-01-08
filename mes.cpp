#include <iostream>

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
const int integralAcc = 1000;       // Dokładność obliczania całki

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
        score += k(x)*u.derivative(x)*v.derivative(x)*h;
    }

    return score;
}

double B(EiFunction u, EiFunction v){
    return calculateIntegral(u, v) - (u.normal(0)*v.normal(0));
}

double L(EiFunction v){
    return -20*v.normal(0);
}


int main(){
    int n;
    cout << "Podaj ilość funkcji bazowych: ";
    cin >> n;
    EiFunction e1(1, n);
    EiFunction e2(2, n);
    cout << B(e1, e1) << endl;
    cout << B(e1, e2) << endl;
    cout << B(e2, e1) << endl;
    cout << B(e2, e2) << endl;

    return 0;
}