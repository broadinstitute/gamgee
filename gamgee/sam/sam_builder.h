#ifndef gamgee__sam_builder__guard
#define gamgee__sam_builder__guard

#include "sam.h"
#include "sam_builder_data_field.h"

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

namespace gamgee {

/**
 * @brief class to build Sam objects from existing data or from scratch
 *
 * Unlike the setters in the Sam class, which only allow in-place modification of
 * the data for efficiency reasons, SamBuilder lets you construct a read in arbitrary
 * ways (including resizing data elements like the number of quality scores and the cigar).
 *
 * To use, either start from scratch with only a header:
 *
 * auto builder = SamBuilder{header};
 *
 * or start with the state in an existing read:
 *
 * auto builder = SamBuilder{existing_read};
 *
 * Then make the desired modifications and build repeatedly.
 *
 * auto read1 = builder.set_cigar("1M1I1M").set_base_quals({1,2,3}).build();
 * auto read2 = builder.set_name("new_name").set_chromosome(3).build();
 *
 * Alternatively, you may use the one_time_build() function if you only need
 * to ever build one read using a given builder. This is slightly more efficient
 * than the more general-purpose repeatable build() function.
 *
 * You must provide a value for all essential data fields before building if you
 * build from scratch, unless you disable validation (not recommended).
 *
 * Some methods of setting values are more efficient than others. For example,
 * setting a cigar using a vector is much more efficient than setting it via a string:
 *
 * builder.set_cigar(vector<CigarElement>{Cigar::make_cigar_element(3, CigarOperator::M)});
 *   is faster than:
 * builder.set_cigar("3M");
 *
 * since the string must be parsed and validated.
 *
 * The builder always copies the data it's given in the set_* methods, even when it belongs
 * to another read, for safety reasons. Altering a read after passing some of its data into
 * a builder will not affect the state of the builder.
 *
 * Before building, the builder performs a validation step to ensure that the read it will
 * create is logically consistent (eg., number of base qualities must match the number of bases).
 * This validation step can be disabled during construction of the builder, but disabling
 * validation is dangerous and not recommended.
 */
class SamBuilder {
 public:
  explicit SamBuilder(const SamHeader& header, const bool validate_on_build = true);                           ///< @brief create a Sam from scratch, starting only with a header
  explicit SamBuilder(const Sam& starting_read, const bool validate_on_build = true);                          ///< @brief create a Sam starting with an existing read and its header
  explicit SamBuilder(const SamHeader& header, const Sam& starting_read, const bool validate_on_build = true); ///< @brief create a Sam starting with an existing read, and manually set the header to a custom value

  // SamBuilders are moveable but not copyable, and use default destruction
  SamBuilder(SamBuilder&& other) = default;
  SamBuilder& operator=(SamBuilder&& other) = default;
  SamBuilder(const SamBuilder& other) = delete;
  SamBuilder& operator=(const SamBuilder& other) = delete;
  ~SamBuilder() = default;

  // Setters for data fields
  SamBuilder& set_name(const std::string& new_name);

  SamBuilder& set_cigar(const Cigar& new_cigar);
  SamBuilder& set_cigar(const std::vector<CigarElement>& new_cigar);
  SamBuilder& set_cigar(const std::initializer_list<CigarElement> new_cigar);
  SamBuilder& set_cigar(const std::string& new_cigar);

  SamBuilder& set_bases(const ReadBases& new_bases);
  SamBuilder& set_bases(const std::vector<Base>& new_bases);
  SamBuilder& set_bases(const std::initializer_list<Base> new_bases);
  SamBuilder& set_bases(const std::string& new_bases);

  SamBuilder& set_base_quals(const BaseQuals& new_base_quals);
  SamBuilder& set_base_quals(const std::vector<uint8_t>& new_base_quals);
  SamBuilder& set_base_quals(const std::initializer_list<uint8_t> new_base_quals);
  SamBuilder& set_base_quals(const std::initializer_list<int> new_base_quals);

  // Setters for core fields: these pass through to the implementations in Sam
  SamBuilder& set_chromosome(const uint32_t chr)              { m_core_read.set_chromosome(chr); return *this;              } ///< @brief simple setter for the chromosome index. Index is 0-based.
  SamBuilder& set_alignment_start(const uint32_t start)       { m_core_read.set_alignment_start(start); return *this;       } ///< @brief simple setter for the alignment start. @warning You should use (1-based and inclusive) alignment but internally this is stored 0-based to simplify BAM conversion.
  SamBuilder& set_mapping_qual(const uint32_t qual)        { m_core_read.set_mapping_qual(qual); return *this;        } ///< @brief simple setter for the mapping quality of the read.
  SamBuilder& set_insert_size(const uint32_t isize)        { m_core_read.set_insert_size(isize); return *this;        } ///< @brief simple setter for the mapping quality of the read.
  SamBuilder& set_mate_chromosome(const uint32_t mchr)        { m_core_read.set_mate_chromosome(mchr); return *this;        } ///< @brief simple setter for the mate's chromosome index. Index is 0-based.
  SamBuilder& set_mate_alignment_start(const uint32_t mstart) { m_core_read.set_mate_alignment_start(mstart); return *this; } ///< @brief simple setter for the mate's alignment start. @warning You should use (1-based and inclusive) alignment but internally this is stored 0-based to simplify BAM conversion.
  SamBuilder& set_paired()            { m_core_read.set_paired(); return *this;            }
  SamBuilder& set_not_paired()        { m_core_read.set_not_paired(); return *this;        }
  SamBuilder& set_properly_paired()            { m_core_read.set_properly_paired(); return *this;            }
  SamBuilder& set_not_properly_paired()        { m_core_read.set_not_properly_paired(); return *this;        }
  SamBuilder& set_unmapped()          { m_core_read.set_unmapped(); return *this;          }
  SamBuilder& set_not_unmapped()      { m_core_read.set_not_unmapped(); return *this;      }
  SamBuilder& set_mate_unmapped()     { m_core_read.set_mate_unmapped(); return *this;     }
  SamBuilder& set_not_mate_unmapped() { m_core_read.set_not_mate_unmapped(); return *this; }
  SamBuilder& set_reverse()           { m_core_read.set_reverse(); return *this;           }
  SamBuilder& set_not_reverse()       { m_core_read.set_not_reverse(); return *this;       }
  SamBuilder& set_mate_reverse()      { m_core_read.set_mate_reverse(); return *this;      }
  SamBuilder& set_not_mate_reverse()  { m_core_read.set_not_mate_reverse(); return *this;  }
  SamBuilder& set_first()             { m_core_read.set_first(); return *this;             }
  SamBuilder& set_not_first()         { m_core_read.set_not_first(); return *this;         }
  SamBuilder& set_last()              { m_core_read.set_last(); return *this;              }
  SamBuilder& set_not_last()          { m_core_read.set_not_last(); return *this;          }
  SamBuilder& set_secondary()         { m_core_read.set_secondary(); return *this;         }
  SamBuilder& set_not_secondary()     { m_core_read.set_not_secondary(); return *this;     }
  SamBuilder& set_fail()              { m_core_read.set_fail(); return *this;              }
  SamBuilder& set_not_fail()          { m_core_read.set_not_fail(); return *this;          }
  SamBuilder& set_duplicate()         { m_core_read.set_duplicate(); return *this;         }
  SamBuilder& set_not_duplicate()     { m_core_read.set_not_duplicate(); return *this;     }
  SamBuilder& set_supplementary()     { m_core_read.set_supplementary(); return *this;     }
  SamBuilder& set_not_supplementary() { m_core_read.set_not_supplementary(); return *this; }

  SamBuilder& clear_tags();                                ///< @brief clear aux tags.
  SamBuilder& add_sam_tags(uint8_t* buffer, const int& len); ///< @brief add tags from hstlib-encoded record.
  SamBuilder& add_char_tag(const std::string& name, const char& value);          ///< @brief add a char-valued tag.
  SamBuilder& add_integer_tag(const std::string& name, const int64_t& value);    ///< @brief add an integer-valued tag.
  SamBuilder& add_float_tag(const std::string& name, const float& value);        ///< @brief add a float-valued tag.
  SamBuilder& add_double_tag(const std::string& name, const double& value);      ///< @brief add a double-valued tag.
  SamBuilder& add_byte_array_tag(const std::string& name, const std::string& value);    ///< @brief add byte array tag.
  SamBuilder& add_string_tag(const std::string& name, const std::string& value);        ///< @brief add a string-valued tag.
  SamBuilder& add_numeric_array_tag(const std::string& name, const SamNumericArrayTag& value); ///< @brief add a numeric array tag.

  Sam build() const;    ///< build a Sam (can be called repeatedly)
  Sam one_time_build(); ///< build a Sam more efficiently by moving the builder's data out of it and invalidating future builds

 private:
  Sam m_core_read;   ///< Shallow Sam object containing only the core (non-data) parts of the read
  SamBuilderDataField m_name;       ///< htslib encoded name for eventual inclusion in the data field
  SamBuilderDataField m_cigar;      ///< htslib encoded cigar for eventual inclusion in the data field
  SamBuilderDataField m_bases;      ///< htslib encoded bases for eventual inclusion in the data field
  SamBuilderDataField m_base_quals; ///< htslib encoded base qualities for eventual inclusion in the data field
  std::unordered_map<std::string, char> m_char_tags;          ////< aux character-valued tags
  std::unordered_map<std::string, int64_t> m_int_tags;        ////< aux integer-valued tags
  std::unordered_map<std::string, float> m_float_tags;        ////< aux float-valued tags
  std::unordered_map<std::string, double> m_double_tags;      ////< aux double-valued tags
  std::unordered_map<std::string, std::string> m_string_tags; ////< aux string-valued tags
  std::unordered_map<std::string, std::string> m_byte_array_tags;           ////< aux byte array tags
  std::unordered_map<std::string, SamNumericArrayTag> m_numeric_array_tags; ////< aux numeric array tags
  bool m_validate_on_build;         ///< should we validate the state of the Sam record at build time?

  void validate() const;                    ///< @brief performs pre-build validation of the state of the Sam record under construction
  void build_tags_array(SamBuilderDataField& tags) const; ///< @brief helper function that constructs the htslib-encoded tags array
  void build_data_array(bam1_t* sam) const; ///< @brief helper function that constructs the concatenated htslib-encoded data array
};

}


#endif /* gamgee__sam_builder__guard */
