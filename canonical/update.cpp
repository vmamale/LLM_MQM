// ============================================================
// update.cpp
// ============================================================
#include "update.h"

void update(Gauge_Field &U) {

  static int    first_time = 1, accept = 0, no_calls = 0;
  static double S_old = 0.0;
  static Gauge_Field f_U;
  static ofstream f_hmc;

  Gauge_Field p_U, old_U, old_f_U;
  double H_old, K_old, K_new, S_new, H_new, hmc_test;

  no_calls++;

  // Refresh momenta with Gaussian distribution
  p_U = Gauge_Field(2);

//  cout << "-------------*****--------------" << endl;
  K_old = kinetic_energy(p_U);

  if (first_time) {
    f_hmc.open("hmc_test");
    if (f_hmc.bad()) cout << "failed to open hmc_test file\n" << flush;
    S_old = action(U);
    force(U, f_U);
    first_time = 0;
  }

  if (no_calls % 1000 == 0) {
    cout << "acceptance rate "
         << (double)accept / (double)no_calls << "\n" << flush;
    no_calls = 0; accept = 0;
  }

  H_old   = S_old + K_old;
  old_U   = U;
  old_f_U = f_U;

//  cout << "Action at beginning of traj " << S_old << "\n" << flush;
//  cout << "H at beginning of traj "      << H_old << "\n" << flush;

  // MD trajectory
  for (int i = 0; i < STEPS; i++)
    evolve_fields(U, p_U, f_U);

  S_new = action(U);
  K_new = kinetic_energy(p_U);
  H_new = S_new + K_new;

//  cout << "Action at end of traj " << S_new << "\n" << flush;
//  cout << "H at end of traj "      << H_new << "\n" << flush;

  // Metropolis accept/reject
  hmc_test = exp(H_old - H_new);
//  f_hmc << hmc_test << "\n" << flush;
//  cout << "exp(-Delta H)	" << hmc_test << endl;

  if (myrandom() < hmc_test) {
    S_old = S_new;
    accept++;
//    cout << "ACCEPTED with deltaH " << H_new - H_old << endl;
   cout << "ACCEPTED" <<  endl;
  } else {
//    cout << "REJECTED with deltaH " << H_new - H_old << endl;
    cout << "REJECTED" << endl;
    U   = old_U;
    f_U = old_f_U;
  }
}
