#include "variant_header.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"

#include <string>
#include <vector>
#include <utility>

#include <iostream>

namespace gamgee {

using namespace std;

/**
 * @brief simple function to plow through all the hrec's and find the ones that match the type tag
 */
vector<string> find_fields_of_type(bcf_hdr_t* h, int type) {
  auto result = vector<string>{};
  for (auto i = 0; i != h->nhrec; ++i) {
    if (h->hrec[i]->type == type)
      result.emplace_back(*(h->hrec[i]->vals));
  }
  return result;
}

VariantHeader::VariantHeader(const VariantHeader& other) :
  m_header {utils::make_shared_variant_header(utils::variant_header_deep_copy(other.m_header.get()))} 
{}

VariantHeader::VariantHeader(VariantHeader&& other) noexcept :
  m_header {move(other.m_header)}
{}

VariantHeader& VariantHeader::operator=(const VariantHeader& other) {
  if ( &other == this )  
    return *this;
  m_header = utils::make_shared_variant_header(utils::variant_header_deep_copy(other.m_header.get())); ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  return *this;
}

/**
 * other is an r-value reference, so it will disappear into the nether right after the swap
 */
VariantHeader& VariantHeader::operator=(VariantHeader&& other) noexcept {
  if (&other == this)
    return *this;
  m_header = move(other.m_header);
  return *this;
}

/**
 * This implementation is simply transforming the char ** representation of the sample names into a 
 * contiguous vector<string>. As efficient as it gets.
 */
vector<string> VariantHeader::samples() const {
  return utils::hts_string_array_to_vector(m_header->samples, uint32_t(bcf_hdr_nsamples(m_header)));
}

vector<string> VariantHeader::contigs() const {
  return find_fields_of_type(m_header.get(), BCF_HL_CTG);
}

vector<string> VariantHeader::filters() const {
  return find_fields_of_type(m_header.get(), BCF_HL_FLT);
}

vector<string> VariantHeader::info_fields() const {
  return find_fields_of_type(m_header.get(), BCF_HL_INFO);
}

vector<string> VariantHeader::format_fields() const {
  return find_fields_of_type(m_header.get(), BCF_HL_FMT);
}

}
