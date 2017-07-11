#include "sam_builder.h"
#include "cigar.h"

#include "htslib/sam.h"

#include "sam_tag.h"
#include "../utils/hts_memory.h"
#include "../utils/sam_tag_utils.h"
#include "../utils/utils.h"

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>

using namespace std;

namespace gamgee {

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
  m_char_tags {},
  m_int_tags {},
  m_float_tags {},
  m_double_tags {},
  m_string_tags {},
  m_byte_array_tags {},
  m_numeric_array_tags {},
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
  m_validate_on_build { validate_on_build }
{
  clear_tags();
  add_sam_tags(bam_get_aux(starting_read.m_body.get()), bam_get_l_aux(starting_read.m_body.get()));
}

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
  m_validate_on_build { validate_on_build }
{
  clear_tags();
  add_sam_tags(bam_get_aux(starting_read.m_body.get()), bam_get_l_aux(starting_read.m_body.get()));
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
  if (new_cigar.empty()) {
    m_cigar.update(nullptr, new_cigar.size() * sizeof(CigarElement), new_cigar.size());
  } else {
    m_cigar.update(&(new_cigar[0]), new_cigar.size() * sizeof(CigarElement), new_cigar.size());
  }
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
  if ( new_cigar.length() == 0 ) {
    m_cigar.update(nullptr, 0, 0);
    return *this;
  }

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
  for ( auto i = 0u; i < num_cigar_elements; ++i ) 
    encoded_cigar_ptr[i] = Cigar::parse_next_cigar_element(new_cigar_stream);

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
 * @brief clear aux tags
 *
 */
SamBuilder& SamBuilder::clear_tags() {
  m_char_tags.clear();
  m_int_tags.clear();
  m_float_tags.clear();
  m_double_tags.clear();
  m_string_tags.clear();
  m_byte_array_tags.clear();
  m_numeric_array_tags.clear();
  return *this;
}

/**
 * @brief add sam tags from hstlib-encoded array
 *
 * @note input tag is added to the list of previously set tags.
 *
 */
SamBuilder& SamBuilder::add_sam_tags(uint8_t* buffer, const int& len) {
  auto tags = utils::parse_encoded_tags(buffer, len);
  for ( const auto& tag : tags ) {
    const auto type = reinterpret_cast<unsigned char>(*tag.second);
    switch (type) {
      case CHAR_SAM_TAG_TYPE_CODE:
        add_char_tag(tag.first, bam_aux2A(tag.second));
        break;
      case INT8_SAM_TAG_TYPE_CODE: case UINT8_SAM_TAG_TYPE_CODE:
      case INT16_SAM_TAG_TYPE_CODE: case UINT16_SAM_TAG_TYPE_CODE:
      case INT32_SAM_TAG_TYPE_CODE: case UINT32_SAM_TAG_TYPE_CODE:
        add_integer_tag(tag.first, bam_aux2i(tag.second));
        break;
      case FLOAT_SAM_TAG_TYPE_CODE:
        add_float_tag(tag.first, bam_aux2f(tag.second));
        break;
      case DOUBLE_SAM_TAG_TYPE_CODE:
        add_double_tag(tag.first, bam_aux2f(tag.second));
        break;
      case STRING_SAM_TAG_TYPE_CODE:
        add_string_tag(tag.first, std::string{bam_aux2Z(tag.second)});
        break;
      case BYTE_ARRAY_SAM_TAG_TYPE_CODE:
        add_byte_array_tag(tag.first, std::string{bam_aux2Z(tag.second)});
        break;
      case NUMERIC_ARRAY_SAM_TAG_TYPE_CODE:
      {
        const auto array_type = reinterpret_cast<unsigned char>(*(tag.second + 1));
        const auto array_len = bam_auxB_len(tag.second);
        switch (array_type) {
          case INT8_SAM_TAG_TYPE_CODE: case UINT8_SAM_TAG_TYPE_CODE:
          case INT16_SAM_TAG_TYPE_CODE: case UINT16_SAM_TAG_TYPE_CODE:
          case INT32_SAM_TAG_TYPE_CODE: case UINT32_SAM_TAG_TYPE_CODE:
          {
            std::vector<int64_t> array_values(array_len);
            for ( uint32_t i = 0; i < array_len; ++i ) {
              array_values[i] = bam_auxB2i(tag.second, i);
            }
            add_numeric_array_tag(
                tag.first,
                SamNumericArrayTag(utils::numeric_array_tag_type(array_type), array_values));
          } break;
          case FLOAT_SAM_TAG_TYPE_CODE:
          {
            std::vector<float> array_values(array_len);
            for ( uint32_t i = 0; i < array_len; ++i ) {
              array_values[i] = bam_auxB2i(tag.second, i);
            }
            add_numeric_array_tag(
                tag.first, SamNumericArrayTag(SamTagType::FLOATARRAY, array_values));
          } break;
          default:
            throw invalid_argument("Unsupported sam numeric array tag type: " +
                                   std::to_string(array_type));
        }
      } break;
      default:
        throw invalid_argument("Unsupported sam tag type: " +
                               std::to_string(type));
    }
  }
  return *this;
}

/**
 * @brief add char-valued tag
 *
 * @note input tag is added to the list of previously set tags.
 *
 */
SamBuilder& SamBuilder::add_char_tag(const std::string& name, const char& value) {
  m_char_tags[name] = value;
  return *this;
}

/**
 * @brief add integer-valued tag
 *
 * @note input tag is added to the list of previously set tags.
 *
 */
SamBuilder& SamBuilder::add_integer_tag(const std::string& name, const int64_t& value) {
  m_int_tags[name] = value;
  return *this;
}

/**
 * @brief add float-valued tag
 *
 * @note input tag is added to the list of previously set tags.
 *
 */
SamBuilder& SamBuilder::add_float_tag(const std::string& name, const float& value) {
  m_float_tags[name] = value;
  return *this;
}

/**
 * @brief add double-valued tag
 *
 * @note input tag is added to the list of previously set tags.
 *
 */
SamBuilder& SamBuilder::add_double_tag(const std::string& name, const double& value) {
  m_double_tags[name] = value;
  return *this;
}

/**
 * @brief add string-valued tag
 *
 * @note input tag is added to the list of previously set tags.
 *
 */
SamBuilder& SamBuilder::add_string_tag(const std::string& name, const string& value) {
  m_string_tags[name] = value;
  return *this;
}

/**
 * @brief add byte array tag
 *
 * @note input tag is added to the list of previously set tags.
 *
 */
SamBuilder& SamBuilder::add_byte_array_tag(const std::string& name, const string& value) {
  m_byte_array_tags[name] = value;
  return *this;
}

/**
 * @brief add numeric array tag
 *
 * @note input tag is added to the list of previously set tags.
 *
 */
SamBuilder& SamBuilder::add_numeric_array_tag(const std::string& name, const SamNumericArrayTag& value) {
  m_numeric_array_tags[name] = value;
  return *this;
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

  // Move m_core_read itself out to the caller, invalidating this builder
  return move(m_core_read);
}

/**
 * @brief performs pre-build validation of the state of the Sam record under construction
 */
void SamBuilder::validate() const {
  // Make sure required data fields have been set
  if ( m_name.empty() || m_bases.empty() || m_base_quals.empty() )
    throw logic_error("Missing one or more required data fields (name, bases, or base qualities)");

  if (!m_core_read.unmapped()) {
    if ( m_cigar.empty() )
      throw logic_error("Missing cigar field");

    // Make sure the sequence length implied by the cigar matches the actual sequence length
    if ( bam_cigar2qlen(m_cigar.num_elements(), (const uint32_t*)m_cigar.raw_data_ptr()) != static_cast<int32_t>(m_bases.num_elements()) )
      throw logic_error("Cigar operations and number of bases do not match");
  }

  // Make sure the number of base qualities matches the number of bases
  if ( m_base_quals.num_elements() != m_bases.num_elements() )
    throw logic_error("Number of bases and number of base qualities do not match");

  // Make sure tag names are 2-character strings
  for ( const auto& tag : m_char_tags ) {
    if (tag.first.length() != SAM_TAG_NAME_LENGTH)
      throw invalid_argument("Tag name should be a " + std::to_string(SAM_TAG_NAME_LENGTH) + "-character string.");
  }
  for ( const auto& tag : m_int_tags ) {
    if (tag.first.length() != SAM_TAG_NAME_LENGTH)
      throw invalid_argument("Tag name should be a " + std::to_string(SAM_TAG_NAME_LENGTH) + "-character string.");
  }
  for ( const auto& tag : m_float_tags ) {
    if (tag.first.length() != SAM_TAG_NAME_LENGTH)
      throw invalid_argument("Tag name should be a " + std::to_string(SAM_TAG_NAME_LENGTH) + "-character string.");
  }
  for ( const auto& tag : m_double_tags ) {
    if (tag.first.length() != SAM_TAG_NAME_LENGTH)
      throw invalid_argument("Tag name should be a " + std::to_string(SAM_TAG_NAME_LENGTH) + "-character string.");
  }
  for ( const auto& tag : m_string_tags ) {
    if (tag.first.length() != SAM_TAG_NAME_LENGTH)
      throw invalid_argument("Tag name should be a " + std::to_string(SAM_TAG_NAME_LENGTH) + "-character string.");
  }
  for ( const auto& tag : m_byte_array_tags ) {
    if (tag.first.length() != SAM_TAG_NAME_LENGTH)
      throw invalid_argument("Tag name should be a " + std::to_string(SAM_TAG_NAME_LENGTH) + "-character string.");
  }
  for ( const auto& tag : m_numeric_array_tags ) {
    if (tag.first.length() != SAM_TAG_NAME_LENGTH)
      throw invalid_argument("Tag name should be a " + std::to_string(SAM_TAG_NAME_LENGTH) + "-character string.");
  }

  // TODO: add additional validation of the core (non-data) fields
}

/**
 * @brief helper function that constructs the concatenated htslib-encoded data array
 */
void SamBuilder::build_data_array(bam1_t* sam) const {
  // Build tags array
  auto tags = SamBuilderDataField();
  build_tags_array(tags);
  const auto data_array_bytes = m_name.num_bytes() + m_cigar.num_bytes() + m_bases.num_bytes() + m_base_quals.num_bytes() + tags.num_bytes();
  // use malloc() instead of new so that htslib can free this memory
  auto data_array = (uint8_t*)malloc(data_array_bytes);
  auto data_array_ptr = data_array;

  data_array_ptr = m_name.copy_into(data_array_ptr);
  data_array_ptr = m_cigar.copy_into(data_array_ptr);
  data_array_ptr = m_bases.copy_into(data_array_ptr);
  data_array_ptr = m_base_quals.copy_into(data_array_ptr);
  data_array_ptr = tags.copy_into(data_array_ptr);

  sam->data = data_array;
  sam->l_data = data_array_bytes;
  sam->m_data = data_array_bytes;
  sam->core.l_qname = m_name.num_bytes();
  sam->core.l_qseq = m_bases.num_elements();
  sam->core.n_cigar = m_cigar.num_elements();
}

/**
 * @brief updates name and type fileds in a htslib-encoded tag array
 */
static uint8_t* update_tag_name_and_type(const std::string& tag_name,
                                         const char& tag_type,
                                         uint8_t* buffer) {
  *buffer++ = tag_name[0];
  *buffer++ = tag_name[1];
  *buffer++ = tag_type;
  return buffer;
}

/**
 * @brief return the size of numeric array values.
 */
static uint32_t numeric_array_value_size(const SamTagType& type) {
  switch (type) {
    case SamTagType::INTEGER8ARRAY:   return sizeof(int8_t);
    case SamTagType::UINTEGER8ARRAY:  return sizeof(uint8_t);
    case SamTagType::INTEGER16ARRAY:  return sizeof(int16_t);
    case SamTagType::UINTEGER16ARRAY: return sizeof(uint16_t);
    case SamTagType::INTEGER32ARRAY:  return sizeof(int32_t);
    case SamTagType::UINTEGER32ARRAY: return sizeof(uint32_t);
    case SamTagType::FLOATARRAY:      return sizeof(float);
    default: throw invalid_argument("Unsupported sam numeric array tag type.");
  }
}

/**
 * @brief helper function that constructs the htslib-encoded tags array
 */
void SamBuilder::build_tags_array(SamBuilderDataField& tags) const {
  // Get number of tags
  const auto num_tags = m_char_tags.size() + m_int_tags.size() + m_float_tags.size() +
      m_double_tags.size() + m_string_tags.size() +
      m_byte_array_tags.size() + m_numeric_array_tags.size();

  auto max_num_bytes = num_tags * ( SAM_TAG_NAME_LENGTH + SAM_TAG_TYPE_LENGTH );
  max_num_bytes += m_char_tags.size() * sizeof(char);
  max_num_bytes += m_int_tags.size() * sizeof(int32_t);
  max_num_bytes += m_float_tags.size() * sizeof(float);
  max_num_bytes += m_double_tags.size() * sizeof(double);
  for ( const auto& tag: m_string_tags )
    max_num_bytes += tag.second.length() + 1;
  for ( const auto& tag: m_byte_array_tags )
    max_num_bytes += tag.second.length() + 1;
  for ( const auto& tag: m_numeric_array_tags ) {
    auto value_size = numeric_array_value_size(tag.second.type());
    max_num_bytes += (tag.second.size() * value_size) +
        // Also, stores type and size of the array.
        sizeof(char) + sizeof(uint32_t);
  }

  auto encoded_tags = std::unique_ptr<uint8_t[]>{ new uint8_t[max_num_bytes] };
  auto encoded_tags_ptr = (uint8_t*) encoded_tags.get();
  // Integer-valued tags could be saved as 1,2,4-byte values depending on
  // their range. "max_num_bytes" reserves the maximum size of 4-byte, but num_bytes
  // adjust the value with the actual number of bytes they are stored at.
  auto num_bytes = max_num_bytes;

  // Append char tags
  // For consistency, append in alphabetical order.
  std::map<std::string, char> ordered_char_tags(
      m_char_tags.cbegin(), m_char_tags.cend());
  for ( const auto& tag : ordered_char_tags ) {
    encoded_tags_ptr = update_tag_name_and_type(tag.first,
                                                CHAR_SAM_TAG_TYPE_CODE,
                                                encoded_tags_ptr);
    *encoded_tags_ptr++ = tag.second;
  }

  // Append integer-valued tags
  // For consistency, append in alphabetical order.
  std::map<std::string, int64_t> ordered_int_tags(
      m_int_tags.cbegin(), m_int_tags.cend());
  for ( const auto& tag : ordered_int_tags ) {
    if ( tag.second >= 0 ) {
      if ( tag.second <= CHAR_MAX ) {
        encoded_tags_ptr = update_tag_name_and_type(
            tag.first, UINT8_SAM_TAG_TYPE_CODE, encoded_tags_ptr);
        encoded_tags_ptr = utils::to_little_endian(uint8_t(tag.second), encoded_tags_ptr);
        num_bytes -= sizeof(int32_t) - sizeof(uint8_t);
      } else if ( tag.second <= USHRT_MAX ) {
        encoded_tags_ptr = update_tag_name_and_type(
            tag.first, UINT16_SAM_TAG_TYPE_CODE, encoded_tags_ptr);
        encoded_tags_ptr = utils::to_little_endian(uint16_t(tag.second), encoded_tags_ptr);
        num_bytes -= sizeof(int32_t) - sizeof(uint16_t);
      } else if ( tag.second <= UINT_MAX ) {
        encoded_tags_ptr = update_tag_name_and_type(
            tag.first, UINT32_SAM_TAG_TYPE_CODE, encoded_tags_ptr);
        encoded_tags_ptr = utils::to_little_endian(uint32_t(tag.second), encoded_tags_ptr);
        num_bytes -= sizeof(int32_t) - sizeof(uint32_t);
      } else throw invalid_argument("Out of range value for an integer tag.");
    } else {
      if ( tag.second >= SCHAR_MIN ) {
        encoded_tags_ptr = update_tag_name_and_type(
            tag.first, INT8_SAM_TAG_TYPE_CODE, encoded_tags_ptr);
        encoded_tags_ptr = utils::to_little_endian(uint8_t(tag.second), encoded_tags_ptr);
        num_bytes -= sizeof(int32_t) - sizeof(int8_t);
      } else if ( tag.second >= SHRT_MIN ) {
        encoded_tags_ptr = update_tag_name_and_type(
            tag.first, INT16_SAM_TAG_TYPE_CODE, encoded_tags_ptr);
        encoded_tags_ptr = utils::to_little_endian(uint16_t(tag.second), encoded_tags_ptr);
        num_bytes -= sizeof(int32_t) - sizeof(int16_t);
      } else if ( tag.second >= INT_MIN ) {
        encoded_tags_ptr = update_tag_name_and_type(
            tag.first, INT32_SAM_TAG_TYPE_CODE, encoded_tags_ptr);
        encoded_tags_ptr = utils::to_little_endian(uint32_t(tag.second), encoded_tags_ptr);
        num_bytes -= sizeof(int32_t) - sizeof(int32_t);
      } else throw invalid_argument("Out of range value for an integer tag.");
    }
  }

  // Append float-valued tags
  // For consistency, append in alphabetical order.
  std::map<std::string, float> ordered_float_tags(
      m_float_tags.cbegin(), m_float_tags.cend());
  for ( const auto& tag : ordered_float_tags ) {
    encoded_tags_ptr = update_tag_name_and_type(tag.first,
                                                FLOAT_SAM_TAG_TYPE_CODE,
                                                encoded_tags_ptr);
    encoded_tags_ptr = utils::to_little_endian(
        *reinterpret_cast<const uint32_t*>(&tag.second),
        encoded_tags_ptr);
  }

  // Append double-valued tags
  // For consistency, append in alphabetical order.
  std::map<std::string, double> ordered_double_tags(
      m_double_tags.cbegin(), m_double_tags.cend());
  for ( const auto& tag : ordered_double_tags ) {
    encoded_tags_ptr = update_tag_name_and_type(tag.first,
                                                DOUBLE_SAM_TAG_TYPE_CODE,
                                                encoded_tags_ptr);
    encoded_tags_ptr = utils::to_little_endian(
        *reinterpret_cast<const uint64_t*>(&tag.second),
        encoded_tags_ptr);
  }

  // Append string tags
  // For consistency, append in alphabetical order.
  std::map<std::string, string> ordered_string_tags(
      m_string_tags.cbegin(), m_string_tags.cend());
  for ( const auto& tag : ordered_string_tags ) {
    encoded_tags_ptr = update_tag_name_and_type(tag.first,
                                                STRING_SAM_TAG_TYPE_CODE,
                                                encoded_tags_ptr);
    memcpy(encoded_tags_ptr, tag.second.c_str(), tag.second.length());
    // Add the null terminator character.
    encoded_tags_ptr[tag.second.length()] = '\0';
    encoded_tags_ptr += tag.second.length() + 1;
  }

  // Append byte array tags
  // For consistency, append in alphabetical order.
  std::map<std::string, string> ordered_byte_array_tags(
      m_byte_array_tags.cbegin(), m_byte_array_tags.cend());
  for ( const auto& tag : ordered_byte_array_tags ) {
    encoded_tags_ptr = update_tag_name_and_type(tag.first,
                                                BYTE_ARRAY_SAM_TAG_TYPE_CODE,
                                                encoded_tags_ptr);
    memcpy(encoded_tags_ptr, tag.second.c_str(), tag.second.length());
    encoded_tags_ptr[tag.second.length()] = '\0';
    encoded_tags_ptr += tag.second.length() + 1;
  }

  // Append numeric array tags
  // For consistency, append in alphabetical order.
  std::map<std::string, SamNumericArrayTag> ordered_numeric_array_tags(
      m_numeric_array_tags.cbegin(), m_numeric_array_tags.cend());
  for ( const auto& tag : ordered_numeric_array_tags ) {
    uint32_t array_len = tag.second.size();
    encoded_tags_ptr = update_tag_name_and_type(tag.first,
                                                NUMERIC_ARRAY_SAM_TAG_TYPE_CODE,
                                                encoded_tags_ptr);
    *encoded_tags_ptr++ = utils::numeric_array_tag_type(tag.second.type());
    encoded_tags_ptr = utils::to_little_endian(array_len, encoded_tags_ptr);
    switch (tag.second.type()) {
      case SamTagType::INTEGER8ARRAY:
        for ( uint32_t i = 0; i < array_len; ++i ) {
          encoded_tags_ptr = utils::to_little_endian(int8_t(tag.second.int_value(i)),
                                                     encoded_tags_ptr);
        }
        break;
      case SamTagType::UINTEGER8ARRAY:
        for ( uint32_t i = 0; i < array_len; ++i ) {
          encoded_tags_ptr = utils::to_little_endian(uint8_t(tag.second.int_value(i)),
                                                     encoded_tags_ptr);
        }
        break;
      case SamTagType::INTEGER16ARRAY:
        for ( uint32_t i = 0; i < array_len; ++i ) {
          encoded_tags_ptr = utils::to_little_endian(int16_t(tag.second.int_value(i)),
                                                     encoded_tags_ptr);
        }
        break;
      case SamTagType::UINTEGER16ARRAY:
        for ( uint32_t i = 0; i < array_len; ++i ) {
          encoded_tags_ptr = utils::to_little_endian(uint16_t(tag.second.int_value(i)),
                                                     encoded_tags_ptr);
        }
        break;
      case SamTagType::INTEGER32ARRAY:
        for ( uint32_t i = 0; i < array_len; ++i ) {
          encoded_tags_ptr = utils::to_little_endian(int32_t(tag.second.int_value(i)),
                                                     encoded_tags_ptr);
        }
        break;
      case SamTagType::UINTEGER32ARRAY:
        for ( uint32_t i = 0; i < array_len; ++i) {
          encoded_tags_ptr = utils::to_little_endian(uint32_t(tag.second.int_value(i)),
                                                     encoded_tags_ptr);
        }
        break;
      case SamTagType::FLOATARRAY:
        for ( uint32_t i = 0; i < array_len; ++i ) {
          float value = tag.second.float_value(i);
          encoded_tags_ptr = utils::to_little_endian(
              *reinterpret_cast<const uint32_t*>(&value),
              encoded_tags_ptr);
        }
        break;
      default:
        throw invalid_argument("Invalid sam numeric array tag type.");
    }
  }

  tags.update(move(encoded_tags), num_bytes, num_tags);
}

}
