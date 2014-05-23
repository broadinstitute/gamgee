#include "sam_pair_iterator.h"
#include "sam.h"
#include "utils/hts_memory.h"

#include "htslib/sam.h"

#include <iostream>
#include <memory>
#include <queue>

using namespace std;

namespace gamgee {

SamPairIterator::SamPairIterator() :
  m_sam_file_ptr    {nullptr},
  m_sam_header_ptr  {nullptr},
  m_sam_record_ptr1 {nullptr},
  m_sam_record_ptr2 {nullptr}
{}

SamPairIterator::SamPairIterator(samFile * sam_file_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr) : 
  m_sam_file_ptr    {sam_file_ptr},
  m_sam_header_ptr  {sam_header_ptr},
  m_sam_record_ptr1 {utils::make_shared_sam(bam_init1())}, ///< important to initialize the record buffer in the constructor so we can reuse it across the iterator
  m_sam_record_ptr2 {utils::make_shared_sam(bam_init1())}, ///< important to initialize the record buffer in the constructor so we can reuse it across the iterator
  m_sam_records     {fetch_next_pair()} ///< important queue must be initialized *before* we call fetch_next_pair. Order matters
{}

SamPairIterator::SamPairIterator(SamPairIterator&& original) :
  m_sam_file_ptr    {original.m_sam_file_ptr},
  m_sam_header_ptr  {move(original.m_sam_header_ptr)},
  m_sam_record_ptr1 {move(original.m_sam_record_ptr1)},
  m_sam_record_ptr2 {move(original.m_sam_record_ptr2)},
  m_sam_records     {move(original.m_sam_records)}
{
  original.m_sam_file_ptr = nullptr;
}

pair<Sam,Sam> SamPairIterator::operator*() {
  return m_sam_records;
}

pair<Sam,Sam> SamPairIterator::operator++() {
  m_sam_records = fetch_next_pair();
  return m_sam_records;
}

bool SamPairIterator::operator!=(const SamPairIterator& rhs) {
  return m_sam_file_ptr != rhs.m_sam_file_ptr;
}

bool SamPairIterator::read_sam(shared_ptr<bam1_t>& record_ptr) {
  if (sam_read1(m_sam_file_ptr, m_sam_header_ptr.get(), record_ptr.get()) < 0) {
    m_sam_file_ptr = nullptr;
    return false;
  }
  return true;
}

Sam SamPairIterator::make_sam(shared_ptr<bam1_t>& record_ptr) {
  return Sam {m_sam_header_ptr, record_ptr};
}

static bool primary(shared_ptr<bam1_t>& record_ptr) {
  return !(record_ptr->core.flag & BAM_FSECONDARY) && !(record_ptr->core.flag & BAM_FSUPPLEMENTARY);
}

Sam SamPairIterator::next_primary_alignment(shared_ptr<bam1_t>& record_ptr) {
  m_supp_alignments.push(utils::make_shared_sam(utils::sam_deep_copy(record_ptr.get())));
  while (read_sam(record_ptr) && !primary(record_ptr))
    m_supp_alignments.push(utils::make_shared_sam(utils::sam_deep_copy(record_ptr.get())));
  return make_sam(record_ptr);
}

pair<Sam,Sam> SamPairIterator::next_supplementary_alignment() {
  const auto read = Sam{m_sam_header_ptr, m_supp_alignments.front()};
  m_supp_alignments.pop();
  return make_pair(read, Sam{});
}

pair<Sam,Sam> SamPairIterator::fetch_next_pair() {
  if (!m_supp_alignments.empty())                                                     // pending supplementary alignments have priority
    return next_supplementary_alignment();
  if (!read_sam(m_sam_record_ptr1))
    return make_pair(Sam{}, Sam{});                                                   // we have reached the end of file
  const auto read1 = make_sam(m_sam_record_ptr1);
  if (!primary(m_sam_record_ptr1) || !read1.paired() || !read_sam(m_sam_record_ptr2)) // unpaired reads go in immediately and by themselves
    return make_pair(read1, Sam{});
  if (primary(m_sam_record_ptr2))                                                     // proper paired alignments return here
    return make_pair(read1, make_sam(m_sam_record_ptr2));
  return make_pair(read1, next_primary_alignment(m_sam_record_ptr2));                 // still haven't found the second primary alignment so search for it while pushing all the secondary/supplementary alignments to the queue
}

}

