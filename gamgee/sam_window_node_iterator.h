#ifndef __gamgee__sam_window_node_iterator__
#define __gamgee__sam_window_node_iterator__

#include "sam.h"
#include "sam_window_node.h"
#include <memory>

namespace gamgee {

class SamWindowNodeIterator {
  public:
    SamWindowNodeIterator();
    SamWindowNodeIterator(SamWindowNode* sam_window_node_ptr, const std::shared_ptr<bam_hdr_t>& sam_header_ptr);

    bool operator!=(SamWindowNodeIterator& rhs);
    Sam& operator*();
    Sam& operator++();

  private:
    std::shared_ptr<bam_hdr_t> m_sam_header_ptr;
    std::shared_ptr<bam1_t> m_sam_record_ptr;
    SamWindowNode* m_sam_window_node_ptr;
    Sam m_sam_record;
    Sam fetch_next_record();
};

}

#endif
