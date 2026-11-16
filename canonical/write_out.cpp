// ============================================================
// write_out.cpp  
// ============================================================
#include "write_out.h"

void write_out(const Gauge_Field &U) {

  static ofstream f_write;

  f_write.clear();                    // reset error flags before reusing
  f_write.open("config");
  if (!f_write.is_open()) {           // correct check: use is_open(), not bad()
    cout << "error opening config file\n" << flush;
    return;                           // don't write into a failed stream
  }

  // Header: T, DT, NCOLOR
  f_write << T << "\t" << DT << "\t" << NCOLOR << "\n";

  for (int t = 0; t < T; t++)
    for (int mu = 0; mu < D; mu++)
      f_write << U.get(t, mu) << "\n";

  f_write << "\n" << flush;
  f_write.close();
}
