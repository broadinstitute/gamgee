#ifndef gamgee__indexed_sam_iterator__guard
#define gamgee__indexed_sam_iterator__guard

#include <iostream> // TODO: Remove post debugging

#include "sam.h"

#include "htslib/sam.h"

#include <memory>
#include <string>
#include <vector>

namespace gamgee {

class IndexedSamIterator {
  public:
    IndexedSamIterator();

    IndexedSamIterator(samFile* sam_file_ptr, hts_idx_t* sam_index_ptr,
        const std::shared_ptr<bam_hdr_t>& sam_header_ptr, const std::vector<std::string>& interval_list);

    IndexedSamIterator(IndexedSamIterator&) = delete;
    IndexedSamIterator& operator=(IndexedSamIterator&) = delete;

    IndexedSamIterator(IndexedSamIterator&&);
    IndexedSamIterator& operator=(IndexedSamIterator&&);

    ~IndexedSamIterator();

    bool operator!=(const IndexedSamIterator& rhs);
    Sam& operator*();
    Sam& operator++();

  private:
    samFile * m_sam_file_ptr;
    hts_idx_t * m_sam_index_ptr;
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr;
    std::vector<std::string> m_interval_list;
    std::vector<std::string>::iterator m_interval_iterator;
    hts_itr_t * m_sam_itr_ptr;
    std::shared_ptr<bam1_t> m_sam_record_ptr;
    Sam m_sam_record;

    hts_itr_t* fetch_current_interval();
    Sam fetch_next_record();
};

}

#endif
