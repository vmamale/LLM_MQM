// ============================================================
// line.cpp  (Polyakov loop — product of U_0 along time)
// ============================================================
#include "line.h"

void line(const Gauge_Field &U) {
  static int first_time = 1, count = 0;
  static ofstream f_line;
  Complex poly;

  count++;
  if (first_time) {
    f_line.open("lines", ios::app);
    if (f_line.bad()) { cerr << "failed to open lines file\n"; exit(1); }
    first_time = 0;
  }

  // Polyakov loop: ordered product of temporal links U_0(t), t=0..T-1
  Umatrix prod = Umatrix(1);
  for (int t = 0; t < T; t++)
    prod = prod * U.get(t, 0);          // mu=0 is the time direction

  poly = (1.0/NCOLOR) * Tr(prod);
//  f_line << poly << "\n" << flush;
  cout << "POLY	" << poly << endl;
}
