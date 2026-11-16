// ============================================================
// setup.cpp  
// ============================================================
#include "setup.h"

void setup(void) {
  cout << "Computing generators for SU(N)\n" << flush;
  (void)my_gen();

  // If RANK == NCOLOR^2, add the U(1) generator
  if (RANK == NCOLOR * NCOLOR) {
    Lambda[RANK-1] = Umatrix();
    for (int i = 0; i < NCOLOR; i++)
      Lambda[RANK-1].set(i, i, (1.0/sqrt(NCOLOR)) * Complex(0.0, 1.0));
  }
}
