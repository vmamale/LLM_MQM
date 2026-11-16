#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

// Parameters
const int NBINS = 200;
const double PI = 3.14159265358979323846;

// Usage: ./hist [input_file] [output_file]
// Input: each line has Re(P) in column 1 and Im(P) in column 2
int main(int argc, char* argv[]) {

  string inName  = (argc > 1) ? argv[1] : "nt20_a0.025/poly";
  string outName = (argc > 2) ? argv[2] : "nt20_a0.025/hist_poly";

  ifstream inputFile(inName);
  if (!inputFile.is_open()) {
    cerr << "Error: Could not open the file " << inName << endl;
    return 1;
  }

  vector<double> histogram(NBINS, 0.0);
  const double dphi = 2 * PI / NBINS;
  long N = 0;  // number of valid samples

  string line;
  while (getline(inputFile, line)) {
    if (line.empty() || line[0] == '#') continue;  // skip blanks/comments

    istringstream iss(line);
    double re, im;
    if (!(iss >> re >> im)) continue;              // need two columns ($1, $2)

    double phase = atan2(im, re);                  // in [-pi, pi]
    int bin = int((phase + PI) / dphi);
    if (bin == NBINS) bin = NBINS - 1;             // phase == +pi edge case
    if (bin >= 0 && bin < NBINS) {
      histogram[bin] += 1.0;
      N++;
    }
  }
  inputFile.close();

  cout << "Read " << N << " Polyakov loop samples\n";  // verify
  if (N == 0) return 1;

  // Output phase density (normalized so integral over (-pi, pi) = 1)
  ofstream out(outName);
  for (int b = 0; b < NBINS; b++) {
    double x   = -PI + (b + 0.5) * dphi;
    double rho = histogram[b] / (N * dphi);
    out << x << " " << rho << "\n";
  }
  out.close();

  return 0;
}
