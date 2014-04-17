#include "sam_pair_iterator.h"
#include "sam.h"

#include "htslib/sam.h"

#include <iostream>
#include <memory>
#include <queue>

using namespace std;

namespace gamgee {

SamPairIterator::SamPairIterator() :
  sam_file_ptr {nullptr},
  sam_header_ptr {nullptr},
  sam_record_ptr {nullptr}
{}

SamPairIterator::SamPairIterator(samFile * sam_file_ptr_, bam_hdr_t * sam_header_ptr_) : 
  sam_file_ptr {sam_file_ptr_},
  sam_header_ptr {sam_header_ptr_},
  sam_record_ptr {bam_init1()}, ///< important to initialize the record buffer in the constructor so we can reuse it across the iterator
  sam_records {fetch_next_pair()} ///< important queue must be initialized *before* we call fetch_next_pair. Order matters
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

Sam SamPairIterator::next_primary_alignment() {
  supp_alignments.emplace(bam_dup1(sam_record_ptr));
  while (read_sam() && not_primary()) 
    supp_alignments.emplace(bam_dup1(sam_record_ptr));
  return make_sam();
}

pair<Sam,Sam> SamPairIterator::next_supplementary_alignment() {
  const auto read = Sam{sam_header_ptr, supp_alignments.front().get()};
  supp_alignments.pop();
  return make_pair(read, Sam{});
}

bool SamPairIterator::not_primary() const {
  return sam_record_ptr->core.flag & BAM_FSECONDARY || sam_record_ptr->core.flag & BAM_FSUPPLEMENTARY;
}

pair<Sam,Sam> SamPairIterator::fetch_next_pair() {
  if (!supp_alignments.empty())                        // pending supplementary alignments have priority
    return next_supplementary_alignment();
  if (!read_sam())
    return make_pair(Sam{}, Sam{});                    // we have reached the end of file
  const auto read1 = make_sam();
  if (not_primary() || !read1.paired() || !read_sam()) // unpaired reads go in immediately and by themselves
    return make_pair(read1, Sam{});
  if (!not_primary())                                  // proper paired alignments return here
    return make_pair(read1, make_sam());
  return make_pair(read1, next_primary_alignment());   // still haven't found the second primary alignment so search for it while pushing all the secondary/supplementary alignments to the queue
}

}

