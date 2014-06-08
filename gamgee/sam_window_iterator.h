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
        samFile* sam_file_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr,
        const uint32_t window_size, const uint32_t step_size);

    SamWindowIterator(SamWindowIterator&&);

    bool operator!=(const SamWindowIterator& rhs);

    SamWindow& operator*();

    SamWindow& operator++();

    ~SamWindowIterator();

  private:
    bool m_is_done_iterating;
    SamWindow m_sam_window; /// sam window
    uint32_t m_window_size;						 /// window size
    uint32_t m_step_size;						 /// step size
    int32_t  m_current_contig;
    int32_t  m_current_position;
    samFile * m_sam_file_ptr;                    ///< pointer to the sam file
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr; ///< pointer to the sam header
    std::shared_ptr<bam1_t> m_sam_record_ptr;    ///< pointer to the internal structure of the sam record. Useful to only allocate it once.

    void move_to_first_window();
    void move_to_next_window();
    bool is_done_sam_window_iterator();
    void delete_sam_window();
    bool is_start_new_window();
    void retrieve_next_record();
    bool is_keep_appending_records(int32_t new_stop_position);
    void append_record();
    bool is_done_pulling_records();
    bool is_keep_dequeueing_records();
    void dequeue_record();
};

}  // end namespace gamgee

#endif
