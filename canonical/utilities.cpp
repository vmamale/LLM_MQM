// ============================================================
// utilities.cpp
// ============================================================
#include "utilities.h"

// ---- Complex ------------------------------------------------
Complex::Complex() { re = 0.0; im = 0.0; }
Complex::Complex(double x, double y) { re = x; im = y; }
double Complex::real() const { return re; }
double Complex::imag() const { return im; }
double Complex::norm() { return sqrt(re*re + im*im); }
void Complex::print() const { cout << "(" << re << ", " << im << ")"; }

ostream& operator<<(ostream& out, Complex c) {
  out << c.real() << "\t" << c.imag();
  return out;
}
istream& operator>>(istream& in, Complex &c) {
  double x, y;
  in >> x >> y;
  c = Complex(x, y);
  return in;
}

Complex operator/(const Complex &o1, const Complex &o2) {
  double norm = o2.real()*o2.real() + o2.imag()*o2.imag();
  return Complex((o1.real()*o2.real() + o1.imag()*o2.imag()) / norm,
                 (o1.imag()*o2.real() - o1.real()*o2.imag()) / norm);
}

Complex pow(const Complex &o1, const int o2) {
  Complex c(1, 0);
  for (int i = 0; i < o2; i++) c = c * o1;
  return c;
}

// ---- Umatrix ------------------------------------------------
Umatrix::Umatrix() {
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      mat[i][j] = Complex();
}

Umatrix::Umatrix(int k) {
  if (k == 1) {
    for (int i = 0; i < NCOLOR; i++) {
      for (int j = 0; j < NCOLOR; j++) mat[i][j] = Complex();
      mat[i][i] = Complex(1.0, 0.0);
    }
  } else
    cout << "wrong Umatrix constructor\n" << flush;
}

Umatrix::Umatrix(Complex m[NCOLOR][NCOLOR]) {
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      mat[i][j] = m[i][j];
}

Complex Umatrix::get(int i, int j) const { return mat[i][j]; }
void    Umatrix::set(int i, int j, const Complex o) { mat[i][j] = o; }

void Umatrix::print() {
  for (int i = 0; i < NCOLOR; i++) {
    for (int j = 0; j < NCOLOR; j++) { mat[i][j].print(); cout << "\t"; }
    cout << "\n";
  }
}

Umatrix Adj(const Umatrix &u) {
  Umatrix res;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      res.set(i, j, conjug(u.get(j, i)));
  return res;
}

ostream& operator<<(ostream& out, Umatrix s) {
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      out << s.get(i, j) << '\t';
  return out;
}
istream& operator>>(istream& in, Umatrix &s) {
  Complex v[NCOLOR][NCOLOR];
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      in >> v[i][j];
  s = Umatrix(v);
  return in;
}

Umatrix operator*(const Umatrix &o1, const Umatrix &o2) {
  Umatrix r;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++) {
      Complex dum;
      for (int k = 0; k < NCOLOR; k++)
        dum = dum + o1.get(i, k) * o2.get(k, j);
      r.set(i, j, dum);
    }
  return r;
}

Umatrix operator+(const Umatrix &x, const Umatrix &y) {
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      dum.set(i, j, x.get(i, j) + y.get(i, j));
  return dum;
}
Umatrix operator-(const Umatrix &x, const Umatrix &y) {
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      dum.set(i, j, x.get(i, j) - y.get(i, j));
  return dum;
}
Umatrix operator*(const Umatrix &o1, const Complex &o2) {
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      dum.set(i, j, o1.get(i, j) * o2);
  return dum;
}
Umatrix operator*(const Complex &o2, const Umatrix &o1) {
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      dum.set(i, j, o1.get(i, j) * o2);
  return dum;
}
Umatrix operator*(const Umatrix &o1, const double o2) {
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      dum.set(i, j, o1.get(i, j) * o2);
  return dum;
}
Umatrix operator*(const double o2, const Umatrix &o1) {
  Umatrix dum;
  for (int i = 0; i < NCOLOR; i++)
    for (int j = 0; j < NCOLOR; j++)
      dum.set(i, j, o1.get(i, j) * o2);
  return dum;
}

Umatrix comm(const Umatrix &o1, const Umatrix &o2) { return o1*o2 - o2*o1; }

Umatrix exp(const Umatrix &u) {
  Umatrix c, del, prod;
  double fac = 1.0;
  int i = 1;
  prod = Umatrix(1);
  c    = Umatrix(1);
  static int sum = 0, counter = 0;

  do {
    fac  = fac * (double)i;
    prod = prod * u;
    del  = prod * (1.0 / fac);
    c    = c + del;
    i++;
  } while (sqrt(Tr(del * Adj(del)).real()) > GAUGETOL && i < 8);

  sum += i; counter++;
  if (counter == 1000) {
    cout << "mean no. of terms in exp() " << (double)sum/counter << "\n" << flush;
    counter = 0; sum = 0;
  }
  return c;
}

Complex Tr(const Umatrix &o) {
  Complex dum;
  for (int i = 0; i < NCOLOR; i++) dum = dum + o.get(i, i);
  return dum;
}

Umatrix real_gaussian_Umatrix() {
  Umatrix dum;
  for (int a = 0; a < RANK; a++)
    dum = dum + gasdev() * Lambda[a];
  return dum;
}

// ---- Gauge_Field --------------------------------------------
// Storage: link[T][D] — T time sites, D=4 links per site

Gauge_Field::Gauge_Field() {
  for (int t = 0; t < T; t++)
    for (int mu = 0; mu < D; mu++)
      link[t][mu] = Umatrix();
}

Gauge_Field::Gauge_Field(int hot) {
  if (hot == 2) {                       // Gaussian: for momenta
    for (int t = 0; t < T; t++)
      for (int mu = 0; mu < D; mu++)
        link[t][mu] = real_gaussian_Umatrix();
    return;
  }
  if (hot == 1) {                       // hot start: random SU(N)
    for (int t = 0; t < T; t++)
      for (int mu = 0; mu < D; mu++)
        link[t][mu] = exp(real_gaussian_Umatrix());
    return;
  }
  if (hot == 0) {                       // cold start: identity
    for (int t = 0; t < T; t++)
      for (int mu = 0; mu < D; mu++)
        link[t][mu] = Umatrix(1);
    return;
  }
  cout << "error in Gauge_Field constructor\n" << flush;
}

Umatrix Gauge_Field::get(const int t, const int mu) const {
  return link[t][mu];
}

void Gauge_Field::set(const int t, const int mu, const Umatrix &u) {
  link[t][mu] = u;
}

Gauge_Field Adj(const Gauge_Field &U) {
  Gauge_Field W;
  for (int t = 0; t < T; t++)
    for (int mu = 0; mu < D; mu++)
      W.set(t, mu, Adj(U.get(t, mu)));
  return W;
}
