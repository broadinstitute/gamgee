#ifndef gamgee__sam__guard
#define gamgee__sam__guard

#include "sam_header.h"
#include "read_bases.h"
#include "base_quals.h"
#include "cigar.h"
#include "sam_tag.h"

#include "htslib/sam.h"

#include <string>
#include <memory>
#include <unordered_map>

namespace gamgee {

/**
 * @brief Utility class to manipulate a Sam record.
 */
class Sam {
 public:
  /**
   * @brief initializes a null Sam.
   * @note this is only used internally by the iterators
   * @warning if you need to create a Sam from scratch, use the builder instead
   */
  Sam() = default;

  /**
   * @brief creates a sam record that points to htslib memory already allocated
   *
   * @note the resulting Sam shares ownership of the pre-allocated memory via shared_ptr
   *       reference counting
   */
  explicit Sam(const std::shared_ptr<bam_hdr_t>& header, const std::shared_ptr<bam1_t>& body) noexcept; 

  /**
   * @brief creates a deep copy of a sam record
   *
   * @note the copy will have exclusive ownership over the newly-allocated htslib memory
   *       until a data field (cigar, bases, etc.) is accessed, after which it will be
   *       shared via reference counting with the Cigar, etc. objects
   * @note does not perform a deep copy of the sam header; to copy the header,
   *       first get it via the header() function and then copy it via the usual C++
   *       semantics
   */
  Sam(const Sam& other);

  /**
   * @copydoc Sam::Sam(const Sam&)  
   */ 
  Sam& operator=(const Sam& other);

  /**
   * @brief moves Sam and its header accordingly. Shared pointers maintain state to all other associated objects correctly.
   */
  Sam(Sam&& other) = default;

  /**
   * @brief move assignment of a Sam and it's header. Shared pointers maintain state to all other associated objects correctly.
   */
  Sam& operator=(Sam&& other) = default;

  /**
   * @brief the header of the Sam record
   *
   * @return a newly created SamHeader object every time it's called but the htslib memory used by the header is the same (no new allocations).
   */
  SamHeader header() const { 
    return SamHeader{m_header}; 
  }

  /**
   * @brief chromosome index of the read.
   *
   * Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getReferenceIndex()
   *
   * @return the integer representation of the chromosome.
   */
  uint32_t chromosome() const { return uint32_t(m_body->core.tid);     }       

  /**
   * @brief the reference position of the first base in the read
   * @note the internal encoding is 0-based to mimic that of the BAM files.
   * @return a (1-based and inclusive) alignment start position (as you would see in a Sam file).
   */
  uint32_t alignment_start() const { return uint32_t(m_body->core.pos+1);   }       

  /**
   * @brief returns a (1-based and inclusive) alignment stop position. 
   * @note the internal encoding is 0-based to mimic that of the BAM files. 
   * @note htslib's bam_endpos returns the coordinate of the first base AFTER the alignment, 0-based, so that translates into the last base IN the 1-based alignment.
   */
  uint32_t alignment_stop() const { return uint32_t(bam_endpos(m_body.get())); }   

  /**
   * @brief calculates the theoretical alignment start of a read that has soft/hard-clips preceding the alignment
   *
   * For example if the read has an alignment start of 100 but the first 4 bases
   * were clipped (hard or soft clipped) then this method will return 96.
   *
   * @return the alignment start (1-based, inclusive) adjusted for clipped bases.  
   *
   * Invalid to call on an unmapped read.
   * Invalid to call with cigar = null
   */
  uint32_t unclipped_start() const ;

  /**
   * @brief calculates the theoretical alignment stop of a read that has soft/hard-clips preceding the alignment
   *
   * For example if the read has an alignment stop of 100 but the last 4 bases
   * were clipped (hard or soft clipped) then this method will return 104.
   *
   * @return the alignment stop (1-based, inclusive) adjusted for clipped bases.     
   * @warning Invalid to call on an unmapped read.
   * @warning Invalid to call with cigar = null
   */
  uint32_t unclipped_stop() const ;

  /**
   * @brief returns the integer representation of the mate's chromosome. 
   *
   * Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). 
   */
  uint32_t mate_chromosome() const { return uint32_t(m_body->core.mtid);    }       

  /**
   * @brief returns a (1-based and inclusive) mate's alignment start position (as you would see in a Sam file). 
   * @note the internal encoding is 0-based to mimic that of the BAM files.
   */
  uint32_t mate_alignment_start() const { return uint32_t(m_body->core.mpos+1);  }       

  /**
   * @brief returns a (1-based and inclusive) mate's alignment stop position.
   * @note the internal encoding is 0-based to mimic that of the BAM files. 
   * @throw std::invalid_argument if called on a record that doesn't contain the mate cigar ("MC") tag.
   */
  uint32_t mate_alignment_stop() const ;                                                

  /**
   * @brief returns a (1-based and inclusive) mate's alignment stop position. 
   *
   * This overload is for usage when the user checks for the existence of the
   * tag themselves and passes it in to avoid exception throwing. This is provided
   * for performance conscious use of this function. This way you will only create 
   * one SamTag object for the mate cigar tag, instead of potentially two when checking
   * for its availability and then calling this function. For example: 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
   * if (!missing(tag))
   *   cout << record.mate_alignment_stop(tag) << endl;  // this will reuse the tag you have already obtained
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * This is better than the alternative using the other overload where you
   * have to either get the Tag twice or check for the exception thrown:
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
   * if (!missing(tag))
   *   cout << record.mate_alignment_stop() << endl;  // this will obtain a new tag internally (unnecessary)
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * @param mate_cigar_tag the MC tag as obtained via the string_tag("MC") API in Sam. 
   * @warning This overload DOES NOT throw an exception if the mate cigar tag is missing. Instead it returns mate_alignment_start(). Treat it as undefined behavior. 
   * @note the internal encoding is 0-based to mimic that of the BAM files. 
   */
  uint32_t mate_alignment_stop(const SamTag<std::string>& mate_cigar_tag) const ;       

  /**
   * @brief returns a (1-based and inclusive) mate's unclipped alignment start position. 
   * @throw std::invalid_argument if called on a record that doesn't contain the mate cigar ("MC") tag.
   */
  uint32_t mate_unclipped_start() const;                                                 

  /**
   * @brief returns a (1-based and inclusive) mate's unclipped alignment start position.
   *
   * This overload is for usage when the user checks for the existence of the
   * tag themselves and passes it in to avoid exception throwing. This is provided
   * for performance conscious use of this function. This way you will only create 
   * one SamTag object for the mate cigar tag. Instead of potentially two when checking
   * for it's availability and then calling this function. For example: 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
   * if (!missing(tag))
   *   cout << record.mate_unclipped_start(tag) << endl;  // this will reuse the tag you have already obtained
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * This is better than the alternative using the other overload where you
   * have to either get the Tag twice or check for the exception thrown:
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
   * if (!missing(tag))
   *   cout << record.mate_unclipped_start() << endl;  // this will obtain a new tag internally (unnecessary)
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * @param mate_cigar_tag the MC tag as obtained via the string_tag("MC") API in Sam. 
   * @warning This overload DOES NOT throw an exception if the mate cigar tag is missing. Instead it returns mate_alignment_start(). Treat it as undefined behavior. 
   * @note the internal encoding is 0-based to mimic that of the BAM files.
   */
  uint32_t mate_unclipped_start(const SamTag<std::string>& mate_cigar_tag) const;        

  /**
   * @brief returns a (1-based and inclusive) mate's unclipped alignment stop position. @throw std::invalid_argument if called on a record that doesn't contain the mate cigar ("MC") tag.
   */
  uint32_t mate_unclipped_stop() const ;                                                

  /**
   * @brief returns a (1-based and inclusive) mate's unclipped alignment stop position. 
   *
   * This overload is for usage when the user checks for the existence of the
   * tag themselves and passes it in to avoid exception throwing. This is provided
   * for performance conscious use of this function. This way you will only create 
   * one SamTag object for the mate cigar tag. Instead of potentially two when checking
   * for it's availability and then calling this function. For example: 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
   * if (!missing(tag))
   *   cout << record.mate_unclipped_stop(tag) << endl;  // this will reuse the tag you have already obtained
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * This is better than the alternative using the other overload where you
   * have to either get the Tag twice or check for the exception thrown:
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
   * if (!missing(tag))
   *   cout << record.mate_unclipped_stop() << endl;  // this will obtain a new tag internally (unnecessary)
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * @param mate_cigar_tag the MC tag as obtained via the string_tag("MC") API in Sam. 
   * @warning This overload DOES NOT throw an exception if the mate cigar tag is missing. Instead it returns mate_alignment_start(). Treat it as undefined behavior. 
   * @note the internal encoding is 0-based to mimic that of the BAM files.
   */
  uint32_t mate_unclipped_stop(const SamTag<std::string>& mate_cigar_tag) const;        

  /**
   * @brief returns the mapping quality of this alignment
   */
  uint8_t mapping_qual() const { return uint8_t(m_body->core.qual); }

  /**
   * @brief inferred insert size as reported by the aligner
   *
   * This is the signed observed insert size. If all segments are mapped to the same reference, 
   * the unsigned observed template length equals the number of bases from the leftmost mapped 
   * base to the rightmost mapped base. The leftmost segment has a plus sign and the rightmost 
   * has a minus sign. The sign of segments in the middle is undefined. 
   *
   * It is set as 0 for single read or when the information is unavailable.
   * 
   * @return a signed insert size or zero if it can't be inferred.
   */
  int32_t insert_size() const { return m_body->core.isize; }

  // modify non-variable length fields (things outside of the data member)
  void set_chromosome(const uint32_t chr)              { m_body->core.tid  = int32_t(chr);        } ///< @brief simple setter for the chromosome index. Index is 0-based.
  void set_alignment_start(const uint32_t start)       { m_body->core.pos  = int32_t(start-1);    } ///< @brief simple setter for the alignment start. @warning You should use (1-based and inclusive) alignment but internally this is stored 0-based to simplify BAM conversion.
  void set_mate_chromosome(const uint32_t mchr)        { m_body->core.mtid = int32_t(mchr);       } ///< @brief simple setter for the mate's chromosome index. Index is 0-based.
  void set_mate_alignment_start(const uint32_t mstart) { m_body->core.mpos = int32_t(mstart - 1); } ///< @brief simple setter for the mate's alignment start. @warning You should use (1-based and inclusive) alignment but internally this is stored 0-based to simplify BAM conversion.
  void set_mapping_qual(const uint8_t mapq)            { m_body->core.qual = mapq;                } ///< @brief simple setter for the alignment quality
  void set_insert_size(const int32_t isize)        { m_body->core.isize = isize;              } ///< @brief simple setter for the insert size

  // getters for fields inside the data field
  std::string name() const { return std::string{bam_get_qname(m_body.get())}; } ///< @brief returns the read name
  Cigar cigar() const { return Cigar{m_body}; }                                 ///< @brief returns the cigar. @warning the objects returned by this member function will share underlying htslib memory with this object. @warning creates an object but doesn't copy the underlying values.
  ReadBases bases() const { return ReadBases{m_body}; }                         ///< @brief returns the read bases. @warning the objects returned by this member function will share underlying htslib memory with this object. @warning creates an object but doesn't copy the underlying values.
  BaseQuals base_quals() const { return BaseQuals{m_body}; }                    ///< @brief returns the base qualities. @warning the objects returned by this member function will share underlying htslib memory with this object. @warning creates an object but doesn't copy the underlying values.

  // getters for tagged values within the aux part of the data field
  SamTag<char> char_tag(const std::string& tag_name) const;          ///< @brief retrieve a char-valued tag by name. @warning creates an object but doesn't copy the underlying values.
  SamTag<int64_t> integer_tag(const std::string& tag_name) const;    ///< @brief retrieve an integer-valued tag by name. @warning creates an object but doesn't copy the underlying values.
  SamTag<double> double_tag(const std::string& tag_name) const;      ///< @brief retrieve an double/float-valued tag by name. @warning creates an object but doesn't copy the underlying values.
  SamTag<std::string> string_tag(const std::string& tag_name) const; ///< @brief retrieve a string-valued tag by name. @warning creates an object but doesn't copy the underlying values.
  SamTag<std::string> byte_array_tag(const std::string& tag_name) const;     ///< @brief retrieve a byte array tag by name. @warning creates an object but doesn't copy the underlying values.
  SamTag<SamNumericArrayTag> numeric_array_tag(const std::string& tag_name) const; ///< @brief retrieve a numeric array tag by name. @warning creates an object but doesn't copy the underlying values.
  std::unordered_map<std::string, SamTagType> all_tags() const; ///< @brief retrieve all tag names and types available in this record
  std::unordered_map<std::string, SamTagType> all_tag_types() const;      ///< @brief retrieve all tag names and types available in this record

  // getters for flags
  bool paired() const { return m_body->core.flag & BAM_FPAIRED;        }          ///< @brief whether or not this read is paired
  bool properly_paired() const { return m_body->core.flag & BAM_FPROPER_PAIR;   } ///< @brief whether or not this read is properly paired (see definition in BAM spec)
  bool unmapped() const { return m_body->core.flag & BAM_FUNMAP;         }        ///< @brief whether or not this read is unmapped
  bool mate_unmapped() const { return m_body->core.flag & BAM_FMUNMAP;        }   ///< @brief whether or not the mate read is unmapped
  bool reverse() const { return m_body->core.flag & BAM_FREVERSE;       }         ///< @brief whether or not this read is from the reverse strand
  bool mate_reverse() const { return m_body->core.flag & BAM_FMREVERSE;      }    ///< @brief whether or not the mate read is from the reverse strand
  bool first() const { return m_body->core.flag & BAM_FREAD1;         }           ///< @brief whether or not this read is the first read in a pair (or multiple pairs)
  bool last() const { return m_body->core.flag & BAM_FREAD2;         }            ///< @brief whether or not this read is the last read in a pair (or multiple pairs)
  bool secondary() const { return m_body->core.flag & BAM_FSECONDARY;     }       ///< @brief whether or not this read is a secondary alignment (see definition in BAM spec)
  bool fail() const { return m_body->core.flag & BAM_FQCFAIL;        }            ///< @brief whether or not this read is marked as failing vendor (sequencer) quality control
  bool duplicate() const { return m_body->core.flag & BAM_FDUP;           }       ///< @brief whether or not this read is a duplicate
  bool supplementary() const { return m_body->core.flag & BAM_FSUPPLEMENTARY; }   ///< @brief whether or not this read is a supplementary alignment (see definition in the BAM spec)

  // modify flags
  void set_paired()            { m_body->core.flag |= BAM_FPAIRED;         }
  void set_not_paired()        { m_body->core.flag &= ~BAM_FPAIRED;        }
  void set_properly_paired()     { m_body->core.flag |= BAM_FPROPER_PAIR;  }
  void set_not_properly_paired() { m_body->core.flag &= ~BAM_FPROPER_PAIR; }
  void set_unmapped()          { m_body->core.flag |= BAM_FUNMAP;          }
  void set_not_unmapped()      { m_body->core.flag &= ~BAM_FUNMAP;         }
  void set_mate_unmapped()     { m_body->core.flag |= BAM_FMUNMAP;         }
  void set_not_mate_unmapped() { m_body->core.flag &= ~BAM_FMUNMAP;        }
  void set_reverse()           { m_body->core.flag |= BAM_FREVERSE;        }
  void set_not_reverse()       { m_body->core.flag &= ~BAM_FREVERSE;       }
  void set_mate_reverse()      { m_body->core.flag |= BAM_FMREVERSE;       }
  void set_not_mate_reverse()  { m_body->core.flag &= ~BAM_FMREVERSE;      }
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

#endif /* defined(gamgee__sam__guard) */
