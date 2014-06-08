/*
 * sam_window_node.cpp
 */

#include "sam_window_node.h"

SamWindowNode::SamWindowNode() {
  next = nullptr;
  sam_record_ptr = nullptr;
}

SamWindowNode::~SamWindowNode() {
  next = nullptr;
  sam_record_ptr = nullptr;
}
