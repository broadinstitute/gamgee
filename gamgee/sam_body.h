#ifndef __gamgee__sam_body__
#define __gamgee__sam_body__

#include <string>

#include "htslib/sam.h"

namespace gamgee {

class SamBody {
 public:
  explicit SamBody();
  SamBody(bam1_t* body);
  SamBody(const SamBody& other);
  SamBody(SamBody&& other);
  SamBody& operator=(const SamBody& other);
  ~SamBody();

  // getters for non-variable length fields (things outside of the data member)
  uint32_t chromosome()           const { return uint32_t(m_body->core.tid);     }
  uint32_t alignment_start()      const { return uint32_t(m_body->core.pos+1);   }
  uint32_t alignment_stop()       const { return uint32_t(bam_endpos(m_body)+1); }
  uint32_t unclipped_start()      const ;
  uint32_t unclipped_stop()       const ;
  uint32_t mate_chromosome()      const { return uint32_t(m_body->core.mtid);    }
  uint32_t mate_alignment_start() const { return uint32_t(m_body->core.mpos+1);  }
  uint32_t mate_alignment_stop()  const ;                                          // not implemented -- requires new mate cigar tag!
  uint32_t mate_unclipped_start() const { return alignment_start();              } // dummy implementation -- requires new mate cigar tag!
  uint32_t mate_unclipped_stop()  const ;                                          // not implemented -- requires new mate cigar tag!

  // modify non-variable length fields (things outside of the data member)
  void set_chromosome(const uint32_t chr)              { m_body->core.tid  = int32_t(chr);        }
  void set_alignment_start(const uint32_t start)       { m_body->core.pos  = int32_t(start-1);    } // incoming alignment is 1-based, storing 0-based
  void set_mate_chromosome(const uint32_t mchr)        { m_body->core.mtid = int32_t(mchr);       }
  void set_mate_alignment_start(const uint32_t mstart) { m_body->core.mpos = int32_t(mstart - 1); } // incoming alignment is 1-based, storing 0-based

  // getters for fields inside the data field
  std::string name()              const { return std::string{bam_get_qname(m_body)}; }

  // getters for flags 
  bool paired()          const { return m_body->core.flag & BAM_FPAIRED;        }
  bool properly_paired() const { return m_body->core.flag & BAM_FPROPER_PAIR;   }
  bool unmapped()        const { return m_body->core.flag & BAM_FUNMAP;         }
  bool next_unmapped()   const { return m_body->core.flag & BAM_FMUNMAP;        }
  bool reverse()         const { return m_body->core.flag & BAM_FREVERSE;       }
  bool next_reverse()    const { return m_body->core.flag & BAM_FMREVERSE;      }
  bool first()           const { return m_body->core.flag & BAM_FREAD1;         }
  bool last()            const { return m_body->core.flag & BAM_FREAD2;         }
  bool secondary()       const { return m_body->core.flag & BAM_FSECONDARY;     }
  bool fail()            const { return m_body->core.flag & BAM_FQCFAIL;        }
  bool duplicate()       const { return m_body->core.flag & BAM_FDUP;           }
  bool supplementary()   const { return m_body->core.flag & BAM_FSUPPLEMENTARY; }

  // modify flags
  void set_paired()            { m_body->core.flag |= BAM_FPAIRED;         }
  void set_not_paired()        { m_body->core.flag &= ~BAM_FPAIRED;        }
  void set_unmapped()          { m_body->core.flag |= BAM_FUNMAP;          }
  void set_not_unmapped()      { m_body->core.flag &= ~BAM_FUNMAP;         }
  void set_next_unmapped()     { m_body->core.flag |= BAM_FMUNMAP;         }
  void set_not_next_unmapped() { m_body->core.flag &= ~BAM_FMUNMAP;        }
  void set_reverse()           { m_body->core.flag |= BAM_FREVERSE;        }
  void set_not_reverse()       { m_body->core.flag &= ~BAM_FREVERSE;       }
  void set_next_reverse()      { m_body->core.flag |= BAM_FMREVERSE;       }
  void set_not_next_reverse()  { m_body->core.flag &= ~BAM_FMREVERSE;      }
  void set_first()             { m_body->core.flag |= BAM_FREAD1;          }
  void set_not_first()         { m_body->core.flag &= ~BAM_FREAD1;         }
  void set_last()              { m_body->core.flag |= BAM_FREAD2;          }
  void set_not_last()          { m_body->core.flag &= ~BAM_FREAD2;         }
  void set_secondary()         { m_body->core.flag |= BAM_FSECONDARY;      }
  void set_not_secondary()     { m_body->core.flag &= ~BAM_FSECONDARY;     }
  void set_fail()              { m_body->core.flag |= BAM_FQCFAIL;         }
  void set_not_fail()          { m_body->core.flag &= ~BAM_FQCFAIL;        }
  void set_duplicate()         { m_body->core.flag |= BAM_FDUP;            }
  void set_not_duplicate()     { m_body->core.flag &= ~BAM_FDUP;           }
  void set_supplementary()     { m_body->core.flag |= BAM_FSUPPLEMENTARY;  }
  void set_not_supplementary() { m_body->core.flag &= ~BAM_FSUPPLEMENTARY; }

  bool empty() const { return m_body->m_data == 0; }

  void make_internal_copy();
 private:
  bam1_t* m_body;
  bool m_must_destroy_body;

  void copy_internal_record(const bam1_t*);

  friend class SamWriter;
};

}
#endif 
