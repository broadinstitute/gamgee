#ifndef __foghorn__sam_body__
#define __foghorn__sam_body__

#include <string>

#include "htslib/sam.h"

namespace foghorn {

class SamBody {
 public:
  explicit SamBody();
  SamBody(const bam1_t* body);
  SamBody(const SamBody& other);
  SamBody(SamBody&& other);
  SamBody& operator=(const SamBody& other);
  ~SamBody();

  // alignment fields
  std::string name()              const { return std::string{bam_get_qname(m_body)}; }
  uint32_t chromosome()           const { return m_body->core.tid; }
  uint32_t alignment_start()      const { return m_body->core.pos+1; }
  uint32_t alignment_stop()       const { return bam_endpos(m_body)+1; }
  uint32_t unclipped_start()      const ;
  uint32_t unclipped_stop()       const ;
  uint32_t mate_chromosome()      const { return m_body->core.mtid; }
  uint32_t mate_alignment_start() const { return m_body->core.mpos+1; }
  uint32_t mate_alignment_stop()  const ;                             // not implemented -- requires new mate cigar tag!
  uint32_t mate_unclipped_start() const { return alignment_start(); } // dummy implementation -- requires new mate cigar tag!
  uint32_t mate_unclipped_stop()  const ;                             // not implemented -- requires new mate cigar tag!

  // check flags 
  bool paired()          const { return m_body->core.flag & BAM_FPAIRED;        }
  bool properly_paired() const { return m_body->core.flag & BAM_FPROPER_PAIR;   }
  bool unmapped()        const { return m_body->core.flag & BAM_FUNMAP;         }
  bool next_unmapped()   const { return m_body->core.flag & BAM_FMUNMAP;        }
  bool first()           const { return m_body->core.flag & BAM_FREAD1;         }
  bool last()            const { return m_body->core.flag & BAM_FREAD2;         }
  bool secondary()       const { return m_body->core.flag & BAM_FSECONDARY;     }
  bool fail()            const { return m_body->core.flag & BAM_FQCFAIL;        }
  bool duplicate()       const { return m_body->core.flag & BAM_FDUP;           }
  bool supplementary()   const { return m_body->core.flag & BAM_FSUPPLEMENTARY; }

  // modify flags
  void mark_as_duplicate() {m_body->core.flag |= BAM_FDUP;}
  void unmark_as_duplicate() {m_body->core.flag &= ~BAM_FDUP;}

  bool empty() const { return m_body->m_data == 0; }

 private:
  bam1_t* m_body;

  friend class SamWriter;
};

}
#endif 
