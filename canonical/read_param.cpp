// ============================================================
// read_param.cpp
//
// Parameter file format (one line):
//   SWEEPS  THERM  GAP SAVE  DT  READIN  SEED
//
// ============================================================
#include "read_param.h"

void read_param(void) {

  ifstream f_in("input");
  if (!f_in.good()) { cout << "\ncan't open file parameters!\n"; exit(1); }

  // File format: SWEEPS THERM GAP SAVE NUMBEROFSTEPS DT READIN SEED
  f_in >> SWEEPS >> THERM >> GAP >> SAVE >> STEPS >> DT >> READIN >> SEED;



   TRAJECTORY_LENGTH = STEPS * DT;

  cout << "One Dimensional reduction of 4D Yang-Mills\n";
  cout << "--------------------------------\n";
  cout << "Number of colors     " << NCOLOR            << "\n";
  cout << "Temporal extent      " << T                 << "\n";
  cout << "Spatial links/site   " << D - 1             << "\n";
  cout << "Lattice Spacing   " << LATSPACE             << "\n";
  cout << "Physical temperature  " << 1.0/(T * LATSPACE) << "\n";
  cout << "Therm sweeps         " << THERM             << "\n";
  cout << "Measurement sweeps   " << SWEEPS            << "\n";
  cout << "Gap                  " << GAP               << "\n";
  cout << "Save                  " << SAVE              << "\n";
  cout << "Time step DT         " << DT                << "\n";
  cout << "Trajectory length    " << TRAJECTORY_LENGTH << "\n";
  cout << "Number of steps    " << STEPS << "\n";
  cout << "Read initial config  " << READIN            << "\n";
  cout << "--------------------------------\n";

  setup();
}
