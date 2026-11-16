// ============================================================
// evolve_fields.cpp
// ============================================================
#include "evolve_fields.h"

void evolve_fields(Gauge_Field &U, Gauge_Field &p_U, Gauge_Field &f_U) {

  Gauge_Field nf_U;

  // Leapfrog step 1: update gauge links
  for (int t = 0; t < T; t++)
    for (int mu = 0; mu < D; mu++)
      U.set(t, mu,
            exp(DT * p_U.get(t,mu) + 0.5 * DT * DT * f_U.get(t,mu))
            * U.get(t,mu));

  // Recompute force at updated U
  force(U, nf_U);

  // Leapfrog step 2: update momenta
  for (int t = 0; t < T; t++)
    for (int mu = 0; mu < D; mu++)
      p_U.set(t, mu,
              p_U.get(t,mu)
              + 0.5 * DT * (nf_U.get(t,mu) + f_U.get(t,mu)));

  // Store new force for next step
  f_U = nf_U;
}
