#ifndef __gamgee__sam_window__
#define __gamgee__sam_window__

#include "sam.h"
#include "sam_window_node.h"

#include "htslib/sam.h"

#include <memory>

namespace gamgee {

/**
 * @brief Utility class to enable for-each style iteration in the SamReader class
 */
class SamWindow {
  public:

    /**
     * @brief creates an empty iterator (used for the end() method)
     */
    SamWindow();

    /**
     * @brief a SamWindow move constructor guarantees all objects will have the same state.
     */
    SamWindow(SamWindow&&);

    /**
     * @brief inequality operator (needed by for-each loop)
     *
     * @param rhs the other SamWindow to compare to
     *
     * @return whether or not the two iterators are the same (e.g. have the same input stream on the same
     * status)
     */
    bool operator!=(const SamWindow& rhs);

    /**
     * @brief dereference operator (needed by for-each loop)
     *
     * @return a persistent Sam object independent from the iterator (a copy of the iterator's object)
     */
    Sam& operator*();

    /**
     * @brief pre-fetches the next record and tests for end of file
     *
     * @return a reference to the object (it can be const& because this return value should only be used
     *         by the for-each loop to check for the eof)
     */
    Sam& operator++();

    std::shared_ptr<bam1_t> peek_record();
    std::shared_ptr<bam1_t> dequeue_record();
    void enqueue_record(std::shared_ptr<bam1_t> sam_record_ptr);
    bool is_empty();

    /**
     * @brief takes care of all the memory allocations of the htslib sam reader interface
     */
    ~SamWindow();

    int32_t m_size;

  private:
    SamWindowNode* m_head;
    SamWindowNode* m_tail;
};

}  // end namespace gamgee

#endif
