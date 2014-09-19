#ifndef gamgee__reference_map__guard
#define gamgee__reference_map__guard

#include <unordered_map>
#include <vector>
#include <string>

#include "fastq_reader.h"
#include "interval.h"

namespace gamgee {

/**
 * @brief Utility class to create a reference object for all reference operations in Foghorn.
 *
 * It hash table stores chr -> sequence as strings internally.
 */
class ReferenceMap : public std::unordered_map<std::string, std::string> {
 public:

  /**
    * @brief opens a fasta/fastq reference file and creates a ReferenceMap from it
    *
    * @param filename reference genome fasta or fastq file
    */
  ReferenceMap(const std::string& filename);

  /**
    * @brief locates the DNA sequence for a given Interval
    *
    * @return DNA sequence for the requested Intervals
    */
  std::string get_sequence(const Interval& interval,         ///< location in the genome
                            const bool reverse_strand = false ///< which strand, releative to the reference genome, to produce the sequence for
      )const;

 private:
  void load_reference_from_fastq_reader(gamgee::FastqReader& reader);
  void read_fastq(gamgee::FastqReader& reader);
};

}  // namespace gamgee

#endif /* gamgee__reference_map__guard */
