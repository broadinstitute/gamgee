#include "sam_builder.h"
#include "utils/hts_memory.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>

using namespace std;

namespace gamgee {

/**
 * @brief Table used to parse chars representing cigar operations into their htslib encodings
 *
 * @note: This table should eventually be moved to htslib. Currently htslib dynamically allocates
 *        and fills the equivalent of this table on-demand, which makes little sense.
 */
const int8_t SamBuilder::cigar_op_parse_table[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,7,-1,-1,-1,-1,9,-1,2,-1,-1,-1,5,1,-1,-1,-1,0,3,-1,6,-1,-1,4,-1,-1,-1,-1,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

/**
 * @brief create a Sam from scratch, starting only with a header
 *
 * @note users of this constructor will need to fill in sufficient parts of the read
 *       to produce a valid read
 */
SamBuilder::SamBuilder(const SamHeader& header, const bool validate_on_build) :
  m_core_read { header.m_header, utils::make_shared_sam(bam_init1()) },
  m_name {},
  m_cigar {},
  m_bases {},
  m_base_quals {},
  m_tags {},
  m_validate_on_build { validate_on_build }
{}

/**
 * @brief create a Sam starting with an existing read and its header
 *
 * @note the data in the existing read is copied into the builder so that future changes
 *       to the read will not affect the builder
 */
SamBuilder::SamBuilder(const Sam& starting_read, const bool validate_on_build) :
  m_core_read { starting_read.m_header, utils::make_shared_sam(utils::sam_shallow_copy(starting_read.m_body.get())) },
  m_name { bam_get_qname(starting_read.m_body.get()), starting_read.m_body.get()->core.l_qname, 1 },
  m_cigar { bam_get_cigar(starting_read.m_body.get()), uint32_t(starting_read.m_body.get()->core.n_cigar << 2), uint32_t(starting_read.m_body.get()->core.n_cigar) },
  m_bases { bam_get_seq(starting_read.m_body.get()), uint32_t((starting_read.m_body.get()->core.l_qseq + 1) >> 1), uint32_t(starting_read.m_body.get()->core.l_qseq) },
  m_base_quals { bam_get_qual(starting_read.m_body.get()), uint32_t(starting_read.m_body.get()->core.l_qseq), uint32_t(starting_read.m_body.get()->core.l_qseq) },
  m_tags { bam_get_aux(starting_read.m_body.get()), uint32_t(bam_get_l_aux(starting_read.m_body.get())), 0 },
  m_validate_on_build { validate_on_build }
{}

/**
 * @brief create a Sam starting with an existing read, and manually set the header to a custom value
 *
 * @note the data in the existing read is copied into the builder so that future changes
 *       to the read will not affect the builder
 */
SamBuilder::SamBuilder(const SamHeader& header, const Sam& starting_read, const bool validate_on_build) :
  m_core_read { header.m_header, utils::make_shared_sam(utils::sam_shallow_copy(starting_read.m_body.get())) },
  m_name { bam_get_qname(starting_read.m_body.get()), starting_read.m_body.get()->core.l_qname, 1 },
  m_cigar { bam_get_cigar(starting_read.m_body.get()), uint32_t(starting_read.m_body.get()->core.n_cigar << 2), uint32_t(starting_read.m_body.get()->core.n_cigar) },
  m_bases { bam_get_seq(starting_read.m_body.get()), uint32_t((starting_read.m_body.get()->core.l_qseq + 1) >> 1), uint32_t(starting_read.m_body.get()->core.l_qseq) },
  m_base_quals { bam_get_qual(starting_read.m_body.get()), uint32_t(starting_read.m_body.get()->core.l_qseq), uint32_t(starting_read.m_body.get()->core.l_qseq) },
  m_tags { bam_get_aux(starting_read.m_body.get()), uint32_t(bam_get_l_aux(starting_read.m_body.get())), 0 },
  m_validate_on_build { validate_on_build }
{}

/**
 * @brief create a builder by move from an existing builder
 */
SamBuilder::SamBuilder(SamBuilder&& other) noexcept :
  m_core_read { move(other.m_core_read) },
  m_name { move(other.m_name) },
  m_cigar { move(other.m_cigar) },
  m_bases { move(other.m_bases) },
  m_base_quals { move(other.m_base_quals) },
  m_tags { move(other.m_tags) },
  m_validate_on_build { other.m_validate_on_build }
{}

/**
 * @brief move an existing builder's state into this object
 */
SamBuilder& SamBuilder::operator=(SamBuilder&& other) noexcept {
  if ( &other == this )
    return *this;

  m_core_read = move(other.m_core_read);
  m_name = move(other.m_name);
  m_cigar = move(other.m_cigar);
  m_bases = move(other.m_bases);
  m_base_quals = move(other.m_base_quals);
  m_tags = move(other.m_tags);
  m_validate_on_build = other.m_validate_on_build;
  return *this;
}

/**
 * @brief set the read's QNAME to the specified value
 */
SamBuilder& SamBuilder::set_name(const std::string& new_name) {
  m_name.update(new_name.c_str(), new_name.length() + 1, 1);  // include a byte for null terminator
  return *this;
}

/**
 * @brief set the read's cigar to the cigar of an existing read
 *
 * @note performs a copy of the other read's cigar
 */
SamBuilder& SamBuilder::set_cigar(const Cigar& new_cigar) {
  m_cigar.update(new_cigar.m_cigar, new_cigar.m_num_cigar_elements * sizeof(CigarElement), new_cigar.m_num_cigar_elements);
  return *this;
}

/**
 * @brief set the read's cigar using a vector of CigarElements
 *
 * @note CigarElements can be created via Cigar::make_cigar_element()
 * @note copies the contents of the vector
 */
SamBuilder& SamBuilder::set_cigar(const std::vector<CigarElement>& new_cigar) {
  m_cigar.update(&(new_cigar[0]), new_cigar.size() * sizeof(CigarElement), new_cigar.size());
  return *this;
}

/**
 * @brief set the read's cigar using an initializer_list of CigarElements
 *
 * @note CigarElements can be created via Cigar::make_cigar_element()
 * @note initializer_lists must be passed by value as per Stroustrup p. 497,
 *       but this is efficient as internally only pointers are copied
 */
SamBuilder& SamBuilder::set_cigar(const std::initializer_list<CigarElement> new_cigar) {
  return set_cigar(vector<CigarElement>{new_cigar});
}

/**
 * @brief set the read's cigar using a string representation of a cigar (eg., "3M1I3M")
 *
 * @note will throw invalid_argument if there is an error parsing the string
 * @note this is the least efficient way to set a cigar
 */
SamBuilder& SamBuilder::set_cigar(const std::string& new_cigar) {
  if ( new_cigar.length() == 0 )
    throw invalid_argument("Empty cigar passed to set_cigar()");

  // Find the number of operators by counting non-digits
  auto num_cigar_elements = new_cigar.length() - count_if(new_cigar.begin(), new_cigar.end(), ::isdigit);
  if ( num_cigar_elements == 0 )
    throw invalid_argument(string("No operators in cigar: ") + new_cigar);

  auto encoded_cigar = unique_ptr<uint8_t[]>{ new uint8_t[num_cigar_elements * sizeof(CigarElement)] };
  auto encoded_cigar_ptr = (uint32_t*)encoded_cigar.get();

  // cannot use auto here since stringstream is not moveable
  // (yeah, "almost always auto" actually involves a lot of overhead in the form of "auto x = type{}"
  // move assignments)
  stringstream new_cigar_stream{ new_cigar };
  for ( auto i = 0u; i < num_cigar_elements; ++i ) {
    unsigned long element_length;
    unsigned char element_op;

    new_cigar_stream >> element_length;
    if ( new_cigar_stream.fail() )
      throw invalid_argument(string("Error parsing cigar string: ") + new_cigar);

    new_cigar_stream >> element_op;
    if ( new_cigar_stream.fail() || int(element_op) >= 128 )
      throw invalid_argument(string("Error parsing cigar string: ") + new_cigar);

    auto encoded_op = cigar_op_parse_table[int(element_op)];
    if ( encoded_op < 0 )
      throw invalid_argument(string("Unrecognized operator ") + char(element_op) + " in cigar string: " + new_cigar);

    encoded_cigar_ptr[i] = (element_length << BAM_CIGAR_SHIFT) | encoded_op;
  }

  // move existing encoded cigar we've just created into the field to avoid an extra copy
  m_cigar.update(move(encoded_cigar), num_cigar_elements * sizeof(CigarElement), num_cigar_elements);
  return *this;
}

/**
 * @brief set the read's bases to the bases of an existing read
 *
 * @note copies the other read's bases
 */
SamBuilder& SamBuilder::set_bases(const ReadBases& new_bases) {
  m_bases.update(new_bases.m_bases, (new_bases.m_num_bases + 1) >> 1, new_bases.m_num_bases);
  return *this;
}

/**
 * @brief set the read's bases using a vector of Bases
 *
 * @note copies the contents of the vector
 */
SamBuilder& SamBuilder::set_bases(const std::vector<Base>& new_bases) {
  auto encoded_size = (new_bases.size() + 1) >> 1;
  auto encoded_bases = unique_ptr<uint8_t[]>{ new uint8_t[encoded_size] };
  auto encoded_base_ptr = encoded_bases.get();

  // 4 bits per base
  memset(encoded_base_ptr, 0, encoded_size);
  for ( auto i = 0u; i < new_bases.size(); ++i ) {
    encoded_base_ptr[i >> 1] |= static_cast<uint8_t>(new_bases[i]) << ((~i & 1) << 2);
  }

  // move the encoded bases we just created into the field to avoid an extra copy
  m_bases.update(move(encoded_bases), encoded_size, new_bases.size());
  return *this;
}

/**
 * @brief set the read's bases using an initializer_list of Bases
 *
 * @note initializer_lists must be passed by value as per Stroustrup p. 497,
 *       but this is efficient as internally only pointers are copied
 */
SamBuilder& SamBuilder::set_bases(const std::initializer_list<Base> new_bases) {
  return set_bases(vector<Base>{new_bases});
}

/**
 * @brief set the read's bases using a string of base values (eg., "ACGT")
 *
 * @note this is the least efficient way to set the bases
 */
SamBuilder& SamBuilder::set_bases(const std::string& new_bases) {
  auto encoded_size = (new_bases.length() + 1) >> 1;
  auto encoded_bases = unique_ptr<uint8_t[]>{ new uint8_t[encoded_size] };
  auto encoded_base_ptr = encoded_bases.get();

  // 4 bits per base
  memset(encoded_base_ptr, 0, encoded_size);
  for ( auto i = 0u; i < new_bases.length(); ++i ) {
    encoded_base_ptr[i >> 1] |= seq_nt16_table[int(new_bases[i])] << ((~i & 1) << 2);
  }

  // move the encoded bases we just created into the field to avoid an extra copy
  m_bases.update(move(encoded_bases), encoded_size, new_bases.length());
  return *this;
}

/**
 * @brief set the read's base qualities to the base qualities of an existing read
 *
 * @note copies the existing read's base qualities
 */
SamBuilder& SamBuilder::set_base_quals(const BaseQuals& new_base_quals) {
  m_base_quals.update(new_base_quals.m_quals, new_base_quals.m_num_quals, new_base_quals.m_num_quals);
  return *this;
}

/**
 * @brief set the read's base qualities using a vector of quality scores
 *
 * @note copies the contents of the vector
 */
SamBuilder& SamBuilder::set_base_quals(const std::vector<uint8_t>& new_base_quals) {
  m_base_quals.update(&(new_base_quals[0]), new_base_quals.size(), new_base_quals.size());
  return *this;
}

/**
 * @brief set the read's base qualities using an initializer_list of (uint8_t) quality scores
 *
 * @note to use this version, you must explicitly specify that the type of elements in the
 *       initializer_list is uint8_t (eg., initializer_list<uint8_t>{1, 2, 3})
 */
SamBuilder& SamBuilder::set_base_quals(const std::initializer_list<uint8_t> new_base_quals) {
  return set_base_quals(vector<uint8_t>{new_base_quals});
}

/**
 * @brief set the read's base qualities using an initializer_list of (int) quality scores
 *
 * @note this version can be used without specifying the type of the elements in
 *       the initializer_list (eg., set_base_quals({1, 2, 3})), but it is less
 *       efficient than the other versions since it must perform a range check on
 *       each score
 */
SamBuilder& SamBuilder::set_base_quals(const std::initializer_list<int> new_base_quals) {
  // Unfortunately, to allow use of initializer lists of quality scores without casting/explicit typing,
  // we need to provide this version that takes an initializer_list<int> instead of an initializer_list<uint8_t>.
  // This means that we must range check the values in the initializer list and explicitly cast them.
  // However, this is acceptable given that this setter is likely to only be called in tests.
  auto quals_vector = vector<uint8_t>(new_base_quals.size());
  auto idx = 0;
  for ( auto qual : new_base_quals ) {
    if ( qual < 0 || qual > 255 )
      throw invalid_argument(string{"Base quality "} + std::to_string(qual) + " invalid. Base qualities must be between 0-255");

    quals_vector[idx] = static_cast<uint8_t>(qual);
    ++idx;
  }
  return set_base_quals(quals_vector);
}

/**
 * @brief build a Sam record using the current state of the builder
 *
 * @note this version of build() can be called repeatedly to build multiple Sam objects with
 *       small variations
 */
Sam SamBuilder::build() const {
  // Perform validation before building, if requested
  if ( m_validate_on_build )
    validate();

  // Allocate memory for a new Sam
  auto new_sam_body = utils::make_shared_sam(bam_init1());
  auto sam_body_ptr = new_sam_body.get();

  // Copy core information first, then construct the concatenated data field
  *(sam_body_ptr) = *(m_core_read.m_body.get());
  build_data_array(sam_body_ptr);

  // Calculate the bin using the same method as htslib
  sam_body_ptr->core.bin = hts_reg2bin(sam_body_ptr->core.pos, sam_body_ptr->core.pos + sam_body_ptr->core.l_qseq, 14, 5);
  // TODO: we don't currently explicitly handle TLEN (core.isize)

  return Sam{ m_core_read.m_header, new_sam_body };
}

/**
 * @brief build a Sam record more efficiently by moving the builder's data out to the caller and invalidating the builder's state
 *
 * @note this method can be called only once, and the builder should not be used after it's called
 */
Sam SamBuilder::one_time_build() {
  // Perform validation before building, if requested
  if ( m_validate_on_build )
    validate();

  // Don't allocate a new Sam -- we'll use the existing memory in m_core_read. We
  // still need to construct the data field, though.
  auto sam_body_ptr = m_core_read.m_body.get();
  build_data_array(sam_body_ptr);

  // Calculate the bin using the same method as htslib
  sam_body_ptr->core.bin = hts_reg2bin(sam_body_ptr->core.pos, sam_body_ptr->core.pos + sam_body_ptr->core.l_qseq, 14, 5);
  // TODO: we don't currently explicitly handle TLEN (core.isize)

  // Move m_core_read itself out to the caller, invalidating this builder
  return move(m_core_read);
}

/**
 * @brief performs pre-build validation of the state of the Sam record under construction
 */
void SamBuilder::validate() const {
  // Make sure required data fields have been set
  if ( m_name.is_empty() || m_cigar.is_empty() || m_bases.is_empty() || m_base_quals.is_empty() )
    throw logic_error("Missing one or more required data fields (name, cigar, bases, or base qualities)");

  // Make sure the sequence length implied by the cigar matches the actual sequence length
  if ( bam_cigar2qlen(m_cigar.num_elements(), (const uint32_t*)m_cigar.raw_data_ptr()) != static_cast<int32_t>(m_bases.num_elements()) )
    throw logic_error("Cigar operations and number of bases do not match");

  // Make sure the number of base qualities matches the number of bases
  if ( m_base_quals.num_elements() != m_bases.num_elements() )
    throw logic_error("Number of bases and number of base qualities do not match");

  // TODO: add additional validation of the core (non-data) fields
}

/**
 * @brief helper function that constructs the concatenated htslib-encoded data array
 */
void SamBuilder::build_data_array(bam1_t* sam) const {
  const auto data_array_bytes = m_name.num_bytes() + m_cigar.num_bytes() + m_bases.num_bytes() + m_base_quals.num_bytes() + m_tags.num_bytes();
  // use malloc() instead of new so that htslib can free this memory
  auto data_array = (uint8_t*)malloc(data_array_bytes);
  auto data_array_ptr = data_array;

  data_array_ptr = m_name.copy_into(data_array_ptr);
  data_array_ptr = m_cigar.copy_into(data_array_ptr);
  data_array_ptr = m_bases.copy_into(data_array_ptr);
  data_array_ptr = m_base_quals.copy_into(data_array_ptr);
  data_array_ptr = m_tags.copy_into(data_array_ptr);

  sam->data = data_array;
  sam->l_data = data_array_bytes;
  sam->m_data = data_array_bytes;
  sam->core.l_qname = m_name.num_bytes();
  sam->core.l_qseq = m_bases.num_elements();
  sam->core.n_cigar = m_cigar.num_elements();
}

}
