#ifndef __gamgee__sam_window_iterator__
#define __gamgee__sam_window_iterator__

#include "sam.h"
#include "sam_window.h"

#include "htslib/sam.h"

#include <memory>

namespace gamgee {

class SamWindowIterator {
  public:

    SamWindowIterator();

    SamWindowIterator(
        samFile* sam_file_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr);

    SamWindowIterator(
        samFile* sam_file_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr,
        const uint32_t window_size, const uint32_t step_size);

    bool operator!=(const SamWindowIterator& rhs);

    SamWindow& operator*();

    SamWindow& operator++();

    bool done_sam_window_iterator() const;

  private:
    samFile * m_sam_file_ptr;
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr;
    std::shared_ptr<bam1_t> m_sam_record_ptr;
    uint32_t m_window_size;
    uint32_t m_step_size;
    int32_t  m_current_contig;
    int32_t  m_current_start;
    int32_t  m_current_stop;
    bool m_done_records;
    SamWindow m_sam_window;

    void move_to_first_window();
    void move_to_next_window();
    void retrieve_next_record();
    bool keep_appending_records();
    void append_record();
    bool keep_dequeuing_records();
    void dequeue_record();
};

}

#endif
