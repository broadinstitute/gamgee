#include "sam_iterator.h"
#include "sam.h"

using namespace std;

namespace gamgee {

SamIterator::SamIterator() :
  sam_file_ptr {nullptr},
  sam_header_ptr {nullptr},
  sam_record_ptr {nullptr}
{}

SamIterator::SamIterator(samFile * sam_file_ptr_, bam_hdr_t * sam_header_ptr_) : 
  sam_file_ptr {sam_file_ptr_},
  sam_header_ptr {sam_header_ptr_},
  sam_record_ptr {bam_init1()},     ///< important to initialize the record buffer in the constructor so we can reuse it across the iterator
  sam_record {fetch_next_record()}
{}

SamIterator::SamIterator(SamIterator&& original) :
  sam_file_ptr {original.sam_file_ptr},
  sam_header_ptr {original.sam_header_ptr},
  sam_record_ptr {original.sam_record_ptr},
  sam_record {original.sam_record}
{}

SamIterator::~SamIterator() {
  bam_destroy1(sam_record_ptr);
  sam_file_ptr = nullptr;
  sam_header_ptr = nullptr;
}

Sam SamIterator::operator*() {
  return sam_record;
}

Sam SamIterator::operator++() {
  sam_record = fetch_next_record();
  return sam_record;
}

bool SamIterator::operator!=(const SamIterator& rhs) {
  return sam_file_ptr != rhs.sam_file_ptr;
}

Sam SamIterator::fetch_next_record() {
  if (sam_read1(sam_file_ptr, sam_header_ptr, sam_record_ptr) < 0) {
    sam_file_ptr = nullptr;
    sam_header_ptr = nullptr;  // can only nullify these two because the SamReader is responsible for freeing them. 
    return Sam{};
  }
  return Sam{sam_header_ptr, sam_record_ptr};
}

}

