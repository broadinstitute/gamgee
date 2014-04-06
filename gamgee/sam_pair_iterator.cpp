#include "sam_pair_iterator.h"
#include "sam.h"

#include <iostream>

using namespace std;

namespace gamgee {

SamPairIterator::SamPairIterator() :
  sam_file_ptr {nullptr},
  sam_header_ptr {nullptr},
  sam_record_ptr {nullptr},
  sam_records {}
{}

SamPairIterator::SamPairIterator(samFile * sam_file_ptr_, bam_hdr_t * sam_header_ptr_) : 
  sam_file_ptr {sam_file_ptr_},
  sam_header_ptr {sam_header_ptr_},
  sam_record_ptr {bam_init1()}, ///< important to initialize the record buffer in the constructor so we can reuse it across the iterator
  sam_records {fetch_next_pair()}
{}

SamPairIterator::SamPairIterator(SamPairIterator&& original) :
  sam_file_ptr {original.sam_file_ptr},
  sam_header_ptr {original.sam_header_ptr},
  sam_record_ptr {original.sam_record_ptr},
  sam_records {original.sam_records}
{
  original.sam_file_ptr = nullptr;
  original.sam_header_ptr = nullptr;
  original.sam_record_ptr = nullptr;
}

SamPairIterator::~SamPairIterator() {
  bam_destroy1(sam_record_ptr);
  sam_file_ptr = nullptr;
  sam_header_ptr = nullptr;
}

pair<Sam,Sam> SamPairIterator::operator*() {
  return sam_records;
}

pair<Sam,Sam> SamPairIterator::operator++() {
  sam_records = fetch_next_pair();
  return sam_records;
}

bool SamPairIterator::operator!=(const SamPairIterator& rhs) {
  return sam_file_ptr != rhs.sam_file_ptr;
}

bool SamPairIterator::read_sam() {
  if (sam_read1(sam_file_ptr, sam_header_ptr, sam_record_ptr) < 0) {
    sam_file_ptr = nullptr;
    sam_header_ptr = nullptr;  // can only nullify these two because the SamReader is responsible for freeing them. 
    return false;
  }
  return true;
}

Sam SamPairIterator::make_sam() {
  return Sam {sam_header_ptr, sam_record_ptr};
}

pair<Sam,Sam> SamPairIterator::fetch_next_pair() {
  if (!read_sam()) 
    return make_pair(Sam{}, Sam{});
  const auto read1 = Sam {sam_header_ptr, sam_record_ptr};
  if (!read1.paired() || !read_sam())
    return make_pair(read1, Sam{});
  const auto read2 = Sam {sam_header_ptr, sam_record_ptr};
  return make_pair(read1, Sam{sam_header_ptr, sam_record_ptr});
}

}

