// ============================================================
// kinetic_energy.cpp
// ============================================================
#include "kinetic_energy.h"

double kinetic_energy(const Gauge_Field &p_U) {
  Complex dum;

  for (int t = 0; t < T; t++)
    for (int mu = 0; mu < D; mu++)
      dum = dum + 0.5 * Tr(Adj(p_U.get(t,mu)) * p_U.get(t,mu));

//  cout << "Umom " << dum.real() << endl;
  return dum.real();
}
