#ifndef gamgee__interval__guard
#define gamgee__interval__guard

#include <vector>
#include <string>

namespace gamgee {

/**
 * @brief Utility class to store an genomic location (Interval).
 *
 * An Interval can be of one or many consecutive bases in a chromosome (contig). We represent
 * Intervals with chromosome, start and stop. A one base Interval has start == stop. No Interval can
 * go beyond a chromosome boundary.
 *
 * The Interval can be output to any output stream (e.g. stdout, file, ...) in any of the supported
 * interval formats. 
 */
class Interval {
 public:

  /** @brief supported output formats for an Interval */
  enum class IntervalType {
      GATK,             ///< GATK format -- chr:start-stop
      PICARD,           ///< PICARD format -- chr [tab] start [tab] stop [tab] strand [tab] comments
      BED               ///< BED format -- chr [tab] start [tab] stop
  };

  /** 
   * @brief creates an interval with standard values
   *
   * chromosome, start and stop are all set to 0. Output type is IntervalType::GATK.
   */
  Interval() :
    m_chr {0},
    m_start {0},
    m_stop {0},
    m_output_type {IntervalType::GATK}
  {}

  /** 
   * @brief creates an interval with all parameters 
   *
   * @param chr         chromosome 
   * @param start       start location (inclusive)
   * @param stop        stop location (inclusive)
   * @param output_type desired output format (only used if Interval is sent to output stream)
   */
  explicit Interval(const std::string& chr, const uint32_t start, const uint32_t stop, const IntervalType output_type = IntervalType::GATK) :
    m_chr{chr},
    m_start{start},
    m_stop{stop},
    m_output_type{output_type} 
  {}

  // it is good to be explicit about these defaults.
  Interval(const Interval&) = default;
  Interval(Interval&&) = default;
  Interval& operator=(const Interval&) = default;
  Interval& operator=(Interval&&) = default;

  /**
    * @brief Checks two Intervals for equality
    *
    * @param rhs the other Interval to compare to
    *
    * @return whether or not the two Intervals have the same chromosome, start and stop. Output
    * type is irrelevant.
    */
  bool operator==(const Interval& rhs) const;

  /**
   * @brief Checks two intervals for inequality
   *
   * Implementation is the negation of operator==
   *
   * @return true if the two intervals don't have the same values in any of their member variables
   */
  bool operator!=(const Interval& rhs) const { return !(*this == rhs); }

  /**
    * @brief calculates the number of loci in the Intervals
    *
    * @return stop() - start() + 1
    */
  inline uint32_t size() const { return m_stop - m_start + 1; }

  /**
    * @brief Tiles an Interval with smaller Intervals anchoring (starting) on the left side (start)
    * of the Interval.
    *
    * Warning: by design, the Interval may or may not be covered on the right (stop) side depending
    * on the length of the original Interval and the parameters passed.
    *
    * @return a vector of Intervals covering the original Interval spaced by spacing loci in
    * between them
    */
  std::vector<Interval> tile_left(const uint32_t tile_size, ///< size of the tiles (small Intervals) to be generated
                                  const uint32_t spacing    ///< how many loci to leave in between the tiles (small Intervals)
                                  ) const;

  /**
    * @brief Tiles an Interval with smaller Intervals anchoring (starting) on the right side (stop)
    * of the Interval.
    *
    * Warning: by design, the Interval may or may not be covered on the left (start) side depending
    * on the length of the original Interval and the parameters passed.
    *
    * @return a vector of Intervals covering the original Interval spaced by spacing loci in
    * between them
    */
  std::vector<Interval> tile_right(const uint32_t tile_size, ///< size of the tiles (small Intervals) to be generated
                                   const uint32_t spacing    ///< how many loci to leave in between the tiles (small Intervals)
                                  ) const;

  /** 
   * @brief creates tiling TSCA baits (small intervals with an insert in between) that tile across the entire Interval with overlaps
   *
   * This tiling procedure will guarantee that the entire Interval is covered (tiled) by creating a
   * list of overlapping tiles that can be used as baits if the design of the overlapping baits is
   * split between the forward and reverse strands of the reference genome. The overlapping bait will
   * always be centered in the spacing between the two baits in the opposite strand.
   *
   * @return a vector of baits / tiles / small intervals that fully cover the original Interval with
   * overlaps.
   */
  std::vector<Interval> tsca_tiling(const uint32_t spacing,     ///< number of loci in between the tiles
                                    const uint32_t insert_size, ///< insert size intended to be used with these baits
                                    const uint32_t flanking     ///< region before and after the target to include in the tiling
                                   );

  std::string chr() const                              { return m_chr;                }
  uint32_t start()  const                              { return m_start;              }
  uint32_t stop()   const                              { return m_stop;               }
  IntervalType output_type() const                     { return m_output_type;        }
  void set_chr(const std::string& chr)                 { m_chr = chr;                 }
  void set_start(const uint32_t start)                 { m_start = start;             }
  void set_stop(const uint32_t stop)                   { m_stop = stop;               }
  void set_output_type(const IntervalType output_type) { m_output_type = output_type; }

 private:
  std::string m_chr;           ///< chromosome
  uint32_t m_start;            ///< first genomic location (inclusive) of the Interval
  uint32_t m_stop;             ///< last genomic location (inclusive) of the Interval
  IntervalType m_output_type;  ///< desired output format for the Interval

};

/**
* @brief utility function to read all Intervals from an Intervals file
*
* @param intervals_file the file with the Intervals in one of the supported formats
*
* @return a vector of Interval objects
*/
std::vector<Interval> read_intervals(const std::string& intervals_file);

/**
* @brief utility function to read all Intervals from an Intervals file
*
* @param input an input stream (e.g. std::cin)
*
* @return a vector of Interval objects
*/
std::vector<Interval> read_intervals(std::istream& input);

}  // end of namespace

/**
* @brief outputs the interval in the outputtype format
*/
std::ostream& operator<<(std::ostream&, const gamgee::Interval&);

#endif /* gamgee__interval__guard */

