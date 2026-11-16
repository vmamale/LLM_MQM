#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
using namespace std;

// ---- Lattice geometry ---------------------------------------
const int T           = 20;
const int D           = 4;
const double LATSPACE = 0.010;

// ---- Gauge group --------------------------------------------
const int NCOLOR = 3;
const int RANK   = NCOLOR * NCOLOR - 1;

// ---- Tolerances ---------------------------------------------
const double GAUGETOL = 1e-6;
const double TRACETOL = 1e-8;

// ---- Global variables (non-class types only) ----------------
extern double DT, TRAJECTORY_LENGTH;
extern int    SWEEPS, GAP, SAVE, THERM, READIN, SEED;
extern int    STEPS;

// ============================================================
// Complex number
// ============================================================
class Complex {
  private:
    double re, im;
  public:
    Complex();
    Complex(double, double);
    double real() const;
    double imag() const;
    double norm();
    void print() const;
    friend ostream& operator<<(ostream &, Complex);
    friend istream& operator>>(istream &, Complex &);
};

inline Complex conjug(const Complex &o) {
  return Complex(o.real(), -o.imag());
}
inline Complex operator+(const Complex &a, const Complex &b) {
  return Complex(a.real()+b.real(), a.imag()+b.imag());
}
inline Complex operator-(const Complex &a, const Complex &b) {
  return Complex(a.real()-b.real(), a.imag()-b.imag());
}
inline Complex operator*(const Complex &a, const Complex &b) {
  return Complex(a.real()*b.real() - a.imag()*b.imag(),
                 a.real()*b.imag() + a.imag()*b.real());
}
inline Complex operator*(const Complex &a, const double b) {
  return Complex(a.real()*b, a.imag()*b);
}
inline Complex operator*(const double a, const Complex &b) {
  return Complex(b.real()*a, b.imag()*a);
}

Complex operator/(const Complex &, const Complex &);
Complex pow(const Complex &, const int);

// ============================================================
// NCOLOR x NCOLOR complex matrix
// ============================================================
class Umatrix {
  private:
    Complex mat[NCOLOR][NCOLOR];
  public:
    Umatrix();
    Umatrix(int);
    Umatrix(Complex [NCOLOR][NCOLOR]);
    Complex get(int, int) const;
    void    set(int, int, const Complex);
    void    print();
    friend ostream& operator<<(ostream &, Umatrix);
    friend istream& operator>>(istream &, Umatrix &);
};

Umatrix operator+(const Umatrix &, const Umatrix &);
Umatrix operator-(const Umatrix &, const Umatrix &);
Umatrix operator*(const Umatrix &, const Umatrix &);
Umatrix operator*(const Umatrix &, const Complex &);
Umatrix operator*(const Complex &, const Umatrix &);
Umatrix operator*(const Umatrix &, const double);
Umatrix operator*(const double,    const Umatrix &);
Umatrix comm(const Umatrix &, const Umatrix &);
Umatrix exp(const Umatrix &);
Umatrix Adj(const Umatrix &);
Complex Tr(const Umatrix &);
Umatrix real_gaussian_Umatrix();

// ---- Now safe to declare: Umatrix is fully defined ----------
extern Umatrix Lambda[RANK];

// ============================================================
// Gauge field: T time sites, D=4 links per site
// ============================================================
class Gauge_Field {
  private:
    Umatrix link[T][D];
  public:
    Gauge_Field();
    Gauge_Field(int);
    Umatrix get(const int t, const int mu) const;
    void    set(const int t, const int mu, const Umatrix &u);
};

Gauge_Field Adj(const Gauge_Field &);

// ============================================================
// Random number utilities
// ============================================================
double gasdev();
double myrandom();

#endif // UTILITIES_H
