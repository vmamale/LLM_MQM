// ============================================================
// measure.cpp
// ============================================================
#include "measure.h"

void measure(const Gauge_Field &U) {
  static int first_time = 1;
  static ofstream f_data;

  if (first_time) {
    f_data.open("data", ios::app);
    if (f_data.bad()) cout << "failed to open data file\n" << flush;
    first_time = 0;
  }

  double act_s;
  double plaq_avg;
  obs(U, act_s, plaq_avg);
//  f_data << act_s / (NCOLOR * NCOLOR * T) << "\n" << flush;
  cout << "Action/(Nc^2*Nt)	" << act_s / (NCOLOR * NCOLOR * T) << endl;
  cout << "RE(PLAQ_AVG)	" << plaq_avg << endl;

  line(U);
}
