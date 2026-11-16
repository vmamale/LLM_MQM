// ============================================================
// obs.cpp  —  observables for the temporal plaquette action
//
// S = (N/2a) * sum_{t=0}^{T-1} sum_{mu=1}^{D-1}
//       Tr( I - U_mu(t) U_0(t) U_mu(t+1)^dag U_0(t)^dag )
//
// Observables computed:
//   act_s      : total action S  (extensive)
//   plaq_avg   : mean plaquette  <Re Tr P / N>,  range [0,1]
//                = 1 - act_s * (2a/N) / (T*(D-1)*N)
//                unity for U=I (weak coupling), 0 at strong coupling
// ============================================================
#include "obs.h"

void obs(const Gauge_Field &U, double &act_s, double &plaq_avg) {

  act_s    = 0.0;
  plaq_avg = 0.0;

  Gauge_Field Udag = Adj(U);

  int n_plaq = T * (D - 1);           // total number of plaquettes

  for (int t = 0; t < T; t++) {
    int tp1 = (t + 1) % T;            // periodic boundary
    for (int mu = 1; mu < D; mu++) {  // spatial links only

      // Plaquette matrix: P(t,mu) = U_mu(t) U_0(t) U_mu(t+1)^dag U_0(t)^dag
      Umatrix P = U.get(t, mu)
                * U.get(t, 0)
                * Udag.get(tp1, mu)
                * Udag.get(t, 0);

      double reTrP = Tr(P).real();

      // Action contribution: Tr(I - P) = N - Re Tr(P)
      act_s    += (double)NCOLOR - reTrP;

      // Plaquette average: Re Tr(P) / N,  normalised to [0,1]
      plaq_avg += reTrP / (double)NCOLOR;
    }
  }

  // Overall prefactor N/(2a) for the action
  act_s    *= (double)NCOLOR / (2.0 * LATSPACE);

  // Normalise plaquette average over all plaquettes
  plaq_avg /= (double)n_plaq;
}
