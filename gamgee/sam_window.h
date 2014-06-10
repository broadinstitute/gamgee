#ifndef __gamgee__sam_window__
#define __gamgee__sam_window__

#include "sam.h"
#include "sam_window_node.h"
#include "sam_window_node_iterator.h"

#include "htslib/sam.h"

#include <memory>

namespace gamgee {

class SamWindow {
  public:

    SamWindow();

    SamWindow(const std::shared_ptr<bam_hdr_t>& sam_header_ptr);

    ~SamWindow();

    std::shared_ptr<bam1_t> peek_record();
    std::shared_ptr<bam1_t> dequeue_record();
    void enqueue_record(std::shared_ptr<bam1_t> sam_record_ptr);
    bool is_empty();

    int32_t start = 0;
    int32_t stop = 0;

    int32_t size() const;

    SamWindowNodeIterator begin() const;
    SamWindowNodeIterator end() const;
  private:
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr;
    SamWindowNode* m_head;
    SamWindowNode* m_tail;
    int32_t m_size = 0;

};

}

#endif
