#include <string>

namespace gamgee {
namespace utils {


/**
 * @brief calculates the reverse complement of a sequence
 * @param sequence the sequence to get the reverse complement from
 * @return a new string with the reverse complement of the sequence
 */
std::string reverse_complement(const std::string& sequence);

/**
 * @brief calculates the complement of a sequence
 * @param sequence the sequence to get the complement from
 * @return a new string with the complement of the sequence
 */
std::string complement(const std::string& sequence); 

/**
 * @brief calculates the complement of a sequence in-place
 * @param sequence the sequence to turn into the complement
 * @return the complement of the sequence calculated in-place (modifies the sequence)
 */
std::string complement(std::string& seq);               

/**
 * @brief calculates the complement of a base
 * @param base the base to get the complement from
 * @return a new char with the complement of base
 */
char complement (const char base);


}
}
