// ============================================================
// check.cpp
// ============================================================
#include "check.h"

void check(Gauge_Field &U) {
  double dummy;

  for (int t = 0; t < T; t++) {
    for (int mu = 0; mu < D; mu++) {
      dummy = fabs((1.0/NCOLOR) * Tr(U.get(t,mu) * Adj(U.get(t,mu))).real() - 1.0);
      if (dummy > GAUGETOL) {
        cout << "U field not unitary - 1/NTr(U*Adj(U))-1 is " << dummy << "\n";
        dummy = sqrt((1.0/NCOLOR) * Tr(U.get(t,mu) * Adj(U.get(t,mu))).real());
        U.set(t, mu, U.get(t,mu) * (1.0/dummy));
      }
    }
  }
}
