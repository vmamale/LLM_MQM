// ============================================================
// read_in.cpp
//
// Reads a saved configuration from the file "config".
//
// Config file header format (must match write_out convention):
//   T  DT  NCOLOR
//
// ============================================================
#include "read_in.h"

void read_in(Gauge_Field &U) {

  int    TDUM;
  double  NCOLORDUM, DTDUM;
  Umatrix dummy;
  ifstream f_read;

  f_read.open("config");
  if (f_read.bad()) cout << "error opening config file\n" << flush;

  // Header: T   DT  NCOLOR
  f_read >> TDUM >> DTDUM >> NCOLORDUM;

  if (TDUM != T)
    cout << "wrong T in config — abort\n" << flush;

  cout << "config time step       " << DTDUM       << "\n";
  cout << "config Ncolor          " << NCOLORDUM   << "\n" << flush;

  // Read all link matrices: loop order must match write_out
  for (int t = 0; t < T; t++)
    for (int mu = 0; mu < D; mu++) {
      f_read >> dummy;
      U.set(t, mu, dummy);
    }

  f_read.close();
}
