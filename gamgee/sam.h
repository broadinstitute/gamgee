#ifndef __gamgee__sam__
#define __gamgee__sam__

#include "sam_header.h"
#include "read_bases.h"
#include "base_quals.h"
#include "cigar.h"
#include "sam_tag.h"

#include "htslib/sam.h"

#include <string>
#include <memory>

namespace gamgee {

/**
 * @brief Utility class to manipulate a Sam record.
 */
class Sam {
 public:
  Sam() = default;                                                                                      ///< @brief initializes a null Sam @note this is only used internally by the iterators @warning if you need to create a Variant from scratch, use the builder instead
  explicit Sam(const std::shared_ptr<bam_hdr_t>& header, const std::shared_ptr<bam1_t>& body) noexcept; ///< @brief creates a Sam given htslib objects. @note used by all iterators
  Sam(const Sam& other);                                                                                ///< @brief makes a deep copy of a Sam and it's header. Shared pointers maintain state to all other associated objects correctly.
  Sam(Sam&& other) noexcept;                                                                            ///< @brief moves Sam and it's header accordingly. Shared pointers maintain state to all other associated objects correctly.
  Sam& operator=(const Sam& other);                                                                     ///< @brief deep copy assignment of a Sam and it's header. Shared pointers maintain state to all other associated objects correctly.
  Sam& operator=(Sam&& other) noexcept;                                                                 ///< @brief move assignment of a Sam and it's header. Shared pointers maintain state to all other associated objects correctly.

  SamHeader header() { return SamHeader{m_header}; }

  // getters for non-variable length fields (things outside of the data member)
  uint32_t chromosome()           const { return uint32_t(m_body->core.tid);     }       ///< @brief returns the integer representation of the chromosome. Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getReferenceIndex()
  uint32_t alignment_start()      const { return uint32_t(m_body->core.pos+1);   }       ///< @brief returns a (1-based and inclusive) alignment start position (as you would see in a Sam file). @note the internal encoding is 0-based to mimic that of the BAM files.
  uint32_t alignment_stop()       const { return uint32_t(bam_endpos(m_body.get())+1); } ///< @brief returns a (1-based and inclusive) alignment stop position (as you would see in a Sam file). @note the internal encoding is 0-based to mimic that of the BAM files.
  uint32_t unclipped_start()      const ;
  uint32_t unclipped_stop()       const ;
  uint32_t mate_chromosome()      const { return uint32_t(m_body->core.mtid);    }       ///< @brief returns the integer representation of the mate's chromosome. Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getMateReferenceIndex()
  uint32_t mate_alignment_start() const { return uint32_t(m_body->core.mpos+1);  }       ///< @brief returns a (1-based and inclusive) mate's alignment start position (as you would see in a Sam file). @note the internal encoding is 0-based to mimic that of the BAM files.
  uint32_t mate_alignment_stop()  const ;                                                ///< @brief returns a (1-based and inclusive) mate'salignment stop position (as you would see in a Sam file). @note the internal encoding is 0-based to mimic that of the BAM files. @warning not implemented -- requires new mate cigar tag!
  uint32_t mate_unclipped_start() const { return alignment_start();              }       ///< @warning dummy implementation -- requires new mate cigar tag!
  uint32_t mate_unclipped_stop()  const ;                                                ///< @warning not implemented -- requires new mate cigar tag!

  // modify non-variable length fields (things outside of the data member)
  // TODO: provide setter for TLEN (core.isize)?
  void set_chromosome(const uint32_t chr)              { m_body->core.tid  = int32_t(chr);        } ///< @brief simple setter for the chromosome index. Index is 0-based.
  void set_alignment_start(const uint32_t start)       { m_body->core.pos  = int32_t(start-1);    } ///< @brief simple setter for the alignment start. @warning You should use (1-based and inclusive) alignment but internally this is stored 0-based to simplify BAM conversion.
  void set_mate_chromosome(const uint32_t mchr)        { m_body->core.mtid = int32_t(mchr);       } ///< @brief simple setter for the mate's chromosome index. Index is 0-based.
  void set_mate_alignment_start(const uint32_t mstart) { m_body->core.mpos = int32_t(mstart - 1); } ///< @brief simple setter for the mate's alignment start. @warning You should use (1-based and inclusive) alignment but internally this is stored 0-based to simplify BAM conversion.

  // getters for fields inside the data field
  std::string name()     const { return std::string{bam_get_qname(m_body.get())}; } ///< @brief returns the read name
  Cigar cigar()          const { return Cigar{m_body}; }                            ///< @brief returns the cigar @warning the objects returned by this member function will share underlying htslib memory with this object
  ReadBases bases()      const { return ReadBases{m_body}; }                        ///< @brief returns the read bases @warning the objects returned by this member function will share underlying htslib memory with this object
  BaseQuals base_quals() const { return BaseQuals{m_body}; }                        ///< @brief returns the base qualities @warning the objects returned by this member function will share underlying htslib memory with this object

  // getters for tagged values within the aux part of the data field
  SamTag<int32_t> integer_tag(const std::string& tag_name) const;    ///< @brief retrieve an integer-valued tag by name
  SamTag<double> double_tag(const std::string& tag_name) const;      ///< @brief retrieve an double/float-valued tag by name
  SamTag<char> char_tag(const std::string& tag_name) const;          ///< @brief retrieve a char-valued tag by name
  SamTag<std::string> string_tag(const std::string& tag_name) const; ///< @brief retrieve a string-valued tag by name

  // getters for flags 
  bool paired()          const { return m_body->core.flag & BAM_FPAIRED;        } ///< @brief whether or not this read is paired
  bool properly_paired() const { return m_body->core.flag & BAM_FPROPER_PAIR;   } ///< @brief whether or not this read is properly paired (see definition in BAM spec)
  bool unmapped()        const { return m_body->core.flag & BAM_FUNMAP;         } ///< @brief whether or not this read is unmapped
  bool next_unmapped()   const { return m_body->core.flag & BAM_FMUNMAP;        } ///< @brief whether or not the next read is unmapped
  bool reverse()         const { return m_body->core.flag & BAM_FREVERSE;       } ///< @brief whether or not this read is from the reverse strand
  bool next_reverse()    const { return m_body->core.flag & BAM_FMREVERSE;      } ///< @brief whether or not the next read is from the reverse strand
  bool first()           const { return m_body->core.flag & BAM_FREAD1;         } ///< @brief whether or not this read is the first read in a pair (or multiple pairs)
  bool last()            const { return m_body->core.flag & BAM_FREAD2;         } ///< @brief whether or not this read is the last read in a pair (or multiple pairs)
  bool secondary()       const { return m_body->core.flag & BAM_FSECONDARY;     } ///< @brief whether or not this read is a secondary alignment (see definition in BAM spec)
  bool fail()            const { return m_body->core.flag & BAM_FQCFAIL;        } ///< @brief whether or not this read is marked as failing vendor (sequencer) quality control
  bool duplicate()       const { return m_body->core.flag & BAM_FDUP;           } ///< @brief whether or not this read is a duplicate
  bool supplementary()   const { return m_body->core.flag & BAM_FSUPPLEMENTARY; } ///< @brief whether or not this read is a supplementary alignment (see definition in the BAM spec) 

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

  bool empty() const { return m_body == nullptr; } ///< @brief whether or not this Sam object is empty, meaning that the internal memory has not been initialized (i.e. a Sam object initialized with Sam()).

 private:
  std::shared_ptr<bam_hdr_t> m_header; ///< htslib pointer to the header structure
  std::shared_ptr<bam1_t> m_body;      ///< htslib pointer to the sam body structure

  friend class SamWriter; ///< allows the writer to access the guts of the object
  friend class SamBuilder; ///< builder needs access to the internals in order to build efficiently
};

}  // end of namespace

#endif /* defined(__gamgee__sam__) */
