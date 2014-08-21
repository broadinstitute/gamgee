#ifndef gamgee__fastq__guard
#define gamgee__fastq__guard

#include <string>

namespace gamgee {

/**
 * @brief Utility class to hold one FastA or FastQ record.
 *
 * Will automatically output a FastA or FastQ based on the presence of quality scores.
 */
class Fastq {

 public:
   
  /** @brief creates an empty record */
  Fastq() :
      m_name {}, m_comment {}, m_sequence{}, m_quals{} {}

  /** @brief creates a full object by assigning all fields */
  Fastq(std::string name,      ///< sequence name
        std::string comment,   ///< optional comment
        std::string sequence,  ///< sequence bases
        std::string quals = "" ///< optional quality scores (leave it out for FastA)
        ) :
      m_name {name}, m_comment {comment}, m_sequence{sequence}, m_quals{quals} 
  {}

  Fastq(const Fastq&) = default;
  Fastq& operator=(const Fastq&) = default;
  Fastq(Fastq&&) = default;
  Fastq& operator=(Fastq&&) = default;

  
  /**
    * @brief inequality comparison of all fields in the record
    *
    * @return true if any field differs (string comparison)
    */
  bool operator!=(const Fastq& other) const {
      return !(*this == other); 
  }
  
  /**
    * @brief equality comparison of all fields in the record
    *
    * @return true only if every field is the same (string comparison)
    */
  bool operator==(const Fastq& other) const {
      return m_name     == other.m_name     &&
             m_comment  == other.m_comment  &&
             m_sequence == other.m_sequence &&
             m_quals    == other.m_quals;
  }


  std::string name() const                       { return m_name;         }
  std::string comment() const                    { return m_comment;      }
  std::string sequence() const                   { return m_sequence;     }
  std::string quals() const                      { return m_quals;        } 
  void set_name(const std::string& name)         { m_name = name;         }
  void set_comment(const std::string& comment)   { m_comment = comment;   }
  void set_sequence(const std::string& sequence) { m_sequence = sequence; }
  void set_quals(const std::string& quals)       { m_quals = quals;       }

  void chop(const int nBases); ///< @brief hard clips the first n bases of the read.
  void reverse_complement();   ///< @brief transform the sequence into it's reverse complement.
  bool is_fastq() const;       ///< @brief true if the record has a quals in it's qual field

 private:

  std::string m_name;     ///< sequence name
  std::string m_comment;  ///< optional comment
  std::string m_sequence; ///< sequence bases
  std::string m_quals;    ///< optional quality scores

};

}  // end of namespace

/**
* @brief outputs the fastq record in fastq format.
*
* The output checks whether the record has quality scores. If it does, it outputs a fastq record,
* otherwise it outputs a fasta record.
*/
std::ostream& operator<< (std::ostream& os, const gamgee::Fastq& fq);

#endif // gamgee__fastq__guard 
