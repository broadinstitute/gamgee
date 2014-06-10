/*
 * sam_window_node.h
 */

#ifndef __gamgee__sam_window_node__
#define __gamgee__sam_window_node__

#include "sam.h"
#include "htslib/sam.h"
#include <memory>

class SamWindowNode {
  public:
    SamWindowNode();
    SamWindowNode* next;
    std::shared_ptr<bam1_t> sam_record_ptr;
};

#endif
