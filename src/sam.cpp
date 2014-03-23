#include <iostream>

#include "sam.h"
namespace gamgee {

Sam::Sam(const bam_hdr_t* header, const bam1_t* body) :
  m_header {header},
  m_body {body} 
{}

Sam::Sam(const Sam& other) :
  m_header {other.m_header},
  m_body {other.m_body}
{}

Sam::Sam(Sam&& other) :
  m_header {other.m_header},
  m_body {other.m_body}
{}

Sam& Sam::operator=(const Sam& other) {
  m_header = other.m_header;
  m_body = other.m_body;
  return *this;
}

} // end of namespace 


using namespace std;
ostream& operator<< (std::ostream& os, const gamgee::Sam& sam) {
  os << sam.chromosome() << ":" << sam.alignment_start() << " | " << sam.mate_chromosome() << ":" << sam.mate_alignment_start();
  return os;
}

