// ============================================================
// main.cpp
//
//compiling instruction: g++ -O3 *.cpp -o hmc -llapack -lblas
//input parameters : input
//run : ./hmc > out
// ============================================================
#include "utilities.h"
#include "update.h"
#include "measure.h"
#include "read_param.h"
#include "read_in.h"
#include "write_out.h"
#include "check.h"
#include "myrandom.h"

int    SWEEPS, GAP, SAVE, THERM, READIN, SEED;
double DT, TRAJECTORY_LENGTH;
int    STEPS;
Umatrix Lambda[RANK];
double SMALLCUT, LARGECUT;

int main(void) {
  int sweep;
  Gauge_Field U;

  SEED = -1;
  read_param();

  if (READIN) {			
    read_in(U);
    check(U);
  } else {
    U = Gauge_Field(1);           // hot start(random) 1, cold start(Identity) 0,
  }				  //  Gaussian: for momenta 2	

  // Thermalization
  cout << "Warming up\n" << flush;
  DT = 0.5 * DT;
  cout << "Thermalization DT " << DT << "\n" << flush;
  for (sweep = 1; sweep <= THERM; sweep++) {
    update(U);
    write_out(U);
  }

  // Measurement sweeps
  DT = 2.0 * DT;
  cout << "Measurement DT " << DT << "\n" << flush;
  cout << "Commencing measurement sweeps\n" << flush;
  for (sweep = 1; sweep <= SWEEPS; sweep++) {
    update(U);
    check(U);
//    measure(U);
//    cout << "sweep no. " << sweep << "\n" << flush;
    if (sweep % GAP == 0) {
      measure(U);
    }
    
    if (sweep % SAVE == 0) {
//      measure(U);
      write_out(U);
    }

  }

  write_out(U);
  return 0;
}
