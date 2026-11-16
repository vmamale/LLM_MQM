// ============================================================
// force.cpp  — HMC forces for the temporal plaquette action
//
// S = (N/2a) * sum_{t,mu>=1} Tr( I - U_mu(t) U_0(t) U_mu(t+1)^dag U_0(t)^dag )
//
// Forces are the left-Lie-algebra derivatives projected to the
// traceless anti-Hermitian (su(N)) part:
//
//   F_U = [Sigma]_ah,    [X]_ah = X - X^dag - (1/N) Tr(X-X^dag) * I
//
// ---- Spatial link U_mu(t), mu >= 1 ----
//
//   U_mu(t) enters P(t,   mu) as the leftmost factor.
//   U_mu(t) enters P(t-1, mu) as U_mu(t)^dag (third factor).
//
//   Sigma_mu(t) = (N/2a) * [
//       - P(t,mu)                                              // forward plaquette
//       + U_mu(t)^dag * U_0(t-1)^dag * U_mu(t-1) * U_0(t-1)  // backward plaquette
//   ]
//
// ---- Temporal link U_0(t) ----
//
//   U_0(t) appears TWICE in each plaquette P(t,mu):
//     (i)  as U_0(t)     [2nd factor]: Tr( U_mu(t) * U_0(t) * B )
//          Lie-deriv of Tr(A U_0 B)  =>  U_0 * B * A
//          term (i)  = U_0(t) * U_mu(t+1)^dag * U_0(t)^dag * U_mu(t)
//
//     (ii) as U_0(t)^dag [4th factor]: Tr( C * U_0(t)^dag )
//          Lie-deriv of Tr(C U^dag)  => -U^dag * C
//          term (ii) = -U_0(t)^dag * U_mu(t) * U_0(t) * U_mu(t+1)^dag
//
//   Both carry an overall minus from -Tr(P), so:
//
//   Sigma_0(t) = -(N/2a) * sum_{mu=1}^{D-1} [
//       U_0(t)     * U_mu(t+1)^dag * U_0(t)^dag * U_mu(t)     <- from U_0
//     - U_0(t)^dag * U_mu(t)       * U_0(t)     * U_mu(t+1)^dag  <- from U_0^dag
//   ]
//
//   Note: after proj_ah (which computes X - X^dag) these two terms reinforce
//   each other, but both must be present before projection.
//
// Both forces are then projected with []_ah and stored in f_U.
// ============================================================
#include "force.h"

// Helper: project to traceless anti-Hermitian (su(N))
static Umatrix proj_ah(const Umatrix &X) {
  Umatrix Y = X - Adj(X);              // anti-Hermitian part
  Complex tr = Tr(Y);
  if (tr.norm() > TRACETOL)
    Y = Y - (1.0 / NCOLOR) * tr * Umatrix(1);
  return Y;
}

void force(const Gauge_Field &U, Gauge_Field &f_U) {

  Gauge_Field Udag = Adj(U);

  for (int t = 0; t < T; t++) {
    int tp1 = (t + 1) % T;
    int tm1 = (t - 1 + T) % T;

    // ----------------------------------------------------------
    // Force on temporal link U_0(t)
    //   U_0(t) appears as both U_0 and U_0^dag in each plaquette,
    //   giving two distinct Lie-derivative contributions per mu:
    //
    //   term_a =  U_0(t)     * U_mu(t+1)^dag * U_0(t)^dag * U_mu(t)
    //   term_b = -U_0(t)^dag * U_mu(t)       * U_0(t)     * U_mu(t+1)^dag
    //
    //   Sigma_0(t) = -(N/2a) * sum_{mu=1}^{D-1} (term_a + term_b)
    // ----------------------------------------------------------
    {
      Umatrix sigma0;
      for (int mu = 1; mu < D; mu++) {
        // term (i): from U_0(t) as the 2nd factor
        Umatrix term_a = U.get(t, 0)
                       * Udag.get(tp1, mu)
                       * Udag.get(t, 0)
                       * U.get(t, mu);

        // term (ii): from U_0(t)^dag as the 4th factor
        Umatrix term_b = Udag.get(t, 0)
                       * U.get(t, mu)
                       * U.get(t, 0)
                       * Udag.get(tp1, mu);
        term_b = term_b * (-1.0);

        sigma0 = sigma0 + term_a + term_b;
      }
      sigma0 = sigma0 * (-(double)NCOLOR / (2.0 * LATSPACE));
      f_U.set(t, 0, proj_ah(sigma0));
    }

    // ----------------------------------------------------------
    // Force on spatial links U_mu(t), mu = 1..D-1
    //   Sigma_mu(t) = (N/2a) * [
    //     - P(t,mu)
    //     + U_mu(t)^dag * U_0(t-1)^dag * U_mu(t-1) * U_0(t-1)
    //   ]
    // where P(t,mu) = U_mu(t) * U_0(t) * U_mu(t+1)^dag * U_0(t)^dag
    // ----------------------------------------------------------
    for (int mu = 1; mu < D; mu++) {
      // Forward plaquette (enters with a minus sign)
      Umatrix P_fwd = U.get(t, mu)
                    * U.get(t, 0)
                    * Udag.get(tp1, mu)
                    * Udag.get(t, 0);

      // Backward plaquette contribution (enters with a plus sign)
      // = U_mu(t)^dag * U_0(t-1)^dag * U_mu(t-1) * U_0(t-1)
      Umatrix P_bwd = Udag.get(t, mu)
                    * Udag.get(tm1, 0)
                    * U.get(tm1, mu)
                    * U.get(tm1, 0);

      Umatrix sigma_mu = ((double)NCOLOR / (2.0 * LATSPACE))
                         * (P_bwd - P_fwd);

      f_U.set(t, mu, proj_ah(sigma_mu));
    }
  }
}
