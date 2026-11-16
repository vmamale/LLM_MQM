// ============================================================
// action.cpp  
//
// S = (N/2a) * sum_{t=0}^{T-1} sum_{mu=1}^{D-1}
//       Tr( I - U_mu(t) * U_0(t) * U_mu(t+1)^dag * U_0(t)^dag )
//
// Layout convention (same as before):
//   mu = 0        : temporal link  U_0(t)
//   mu = 1..D-1   : spatial links  U_mu(t)
//
// Time is periodic: t+1 is taken mod T.
// ============================================================
#include "action.h"

double action(const Gauge_Field &U) {
  double act_s = 0.0;
  Gauge_Field Udag = Adj(U);

  for (int t = 0; t < T; t++) {
    int tp1 = (t + 1) % T;                   // periodic boundary
    for (int mu = 1; mu < D; mu++) {          // spatial links only

      Umatrix P = U.get(t, mu)
                * U.get(t, 0)
                * Udag.get(tp1, mu)
                * Udag.get(t, 0);

      // Tr(I - P) = N - Re Tr(P)
      act_s += (double)NCOLOR - Tr(P).real();
    }
  }

  return ((double)NCOLOR / (2.0 * LATSPACE)) * act_s;
}
