#ifndef __gamgee__sam__
#define __gamgee__sam__

#include <string>

#include "htslib/sam.h"
#include "sam_header.h"
#include "sam_body.h"

namespace gamgee {

/**
 * @brief Utility class to manipulate a Sam record.
 */
class Sam {
 public:
  explicit Sam() = default;
  Sam(const bam_hdr_t* header, const bam1_t* body);
  Sam(const Sam& other);
  Sam(Sam&& other);
  Sam& operator=(const Sam& other);

  // alignment fields
  std::string name()              const { return m_body.name(); }
  uint32_t chromosome()           const { return m_body.chromosome(); }
  uint32_t alignment_start()      const { return m_body.alignment_start(); }
  uint32_t alignment_stop()       const { return m_body.alignment_stop(); }
  uint32_t unclipped_start()      const { return m_body.unclipped_start(); }
  uint32_t unclipped_stop()       const { return m_body.unclipped_stop(); }
  uint32_t mate_chromosome()      const { return m_body.mate_chromosome(); }
  uint32_t mate_alignment_start() const { return m_body.mate_alignment_start(); }

  // check flags 
  bool paired()          const { return m_body.paired();          }
  bool unpaired()        const { return !m_body.paired();         }
  bool properly_paired() const { return m_body.properly_paired(); }
  bool mapped()          const { return !m_body.unmapped();       }
  bool unmapped()        const { return m_body.unmapped();        }
  bool mate_mapped()     const { return !m_body.next_unmapped();  }
  bool mate_unmapped()   const { return m_body.next_unmapped();   }
  bool first()           const { return m_body.first();           }
  bool last()            const { return m_body.last();            }
  bool secondary()       const { return m_body.secondary();       }
  bool fail()            const { return m_body.fail();            }
  bool duplicate()       const { return m_body.duplicate();       }
  bool unique()          const { return !m_body.duplicate();      }
  bool supplementary()   const { return m_body.supplementary();   }

  // modify flags
  void mark_as_duplicate()   { m_body.mark_as_duplicate();}
  void unmark_as_duplicate() { m_body.unmark_as_duplicate();}

  bool empty() const { return m_body.empty(); }

  /**
   * @brief shares the private member data as const to prohibit others from modifying 
   * while still allowing access. If you need to modify it, make a copy with the provided
   * SamHeader(const SamHeader&) copy constructor
   *
   * @return the header of this record
   */
  const SamHeader& header() const { return m_header; }

  /**
   * @brief shares the private member data as const to prohibit others from modifying 
   * while still allowing access. If you need to modify it, make a copy with the provided
   * SamBody(const SamBody&) copy constructor
   *
   * @return the body of this record
   */
  const SamBody& body() const { return m_body; }


 private:
  SamHeader m_header;
  SamBody m_body;

};

}  // end of namespace

/**
 * @brief outputs the sam record in sam format.
 *
 * The output checks whether the record has quality scores. If it does, it outputs a sam record,
 * otherwise it outputs a fasta record.
 */
std::ostream& operator<< (std::ostream& os, const gamgee::Sam& sam);

#endif /* defined(__gamgee__sam__) */
