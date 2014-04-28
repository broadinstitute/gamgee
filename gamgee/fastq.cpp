#include "fastq.h"
#include "utils.h"

#include <iostream>

using namespace std;

namespace gamgee { 

bool Fastq::is_fastq() const {
  return !m_quals.empty();
}

void Fastq::chop(const int nBases) {
  m_sequence.erase(0, nBases); 
  if (is_fastq()) 
    m_quals.erase(0, nBases);
}

void Fastq::reverse_complement() {
  m_sequence = utils::reverse_complement(m_sequence);
}

}  // end of namespace

static string print_fasta_record(const gamgee::Fastq& fq) {
  return ">" + fq.name() + " " + fq.comment() + "\n" + fq.sequence() + "\n";
}

static string print_fastq_record(const gamgee::Fastq& fq) {
  return "@" + fq.name() + " " + fq.comment() + "\n" + fq.sequence() + "\n+\n" + fq.quals() + "\n";
}

std::ostream& operator<< (std::ostream& os, const gamgee::Fastq& fq) {
  return os << (fq.is_fastq() ? print_fastq_record(fq) : print_fasta_record(fq));
}
