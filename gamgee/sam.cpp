#include "sam.h"

#include <iostream>

using namespace std;
ostream& operator<< (std::ostream& os, const gamgee::Sam& sam) {
  os << sam.chromosome() << ":" << sam.alignment_start() << " | " << sam.mate_chromosome() << ":" << sam.mate_alignment_start();
  return os;
}

