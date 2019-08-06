#include "variant_header.h"

#include "../utils/hts_memory.h"
#include "../utils/utils.h"

#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include <iostream>

namespace gamgee {

using namespace std;

/**
 * @brief simple function to plow through all the hrec's and find the ones that match the type tag
 */
static vector<string> find_fields_of_type(bcf_hdr_t* h, int type) {
  auto result = vector<string>{};
  for (auto i = 0; i != h->nhrec; ++i) {
    if (h->hrec[i]->type == type)
      result.emplace_back(*(h->hrec[i]->vals));
  }
  return result;
}

static uint32_t count_fields_of_type(const bcf_hdr_t* h, const int type) {
  auto count = 0u;
  for ( auto i = 0; i != h->nhrec; ++i ) {
    if ( h->hrec[i]->type == type )
      ++count;
  }
  return count;
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

bool VariantHeader::operator==(const VariantHeader& rhs) const {
  // compare names
  if (samples() != rhs.samples()) return false;
  if (chromosomes() != rhs.chromosomes()) return false;
  if (filters() != rhs.filters()) return false;
  if (shared_fields() != rhs.shared_fields()) return false;
  if (individual_fields() != rhs.individual_fields()) return false;

  // can't use index here because filters/shared/individual indices depend on the insertion order of the others

  for (const auto& field : shared_fields()) {
    if (shared_field_type(field) != rhs.shared_field_type(field))
      return false;
  }

  for (const auto& field : individual_fields()) {
    if (individual_field_type(field) != rhs.individual_field_type(field))
      return false;
  }

  // TODO? This covers everything that VariantHeader provides access to,
  // but the underlying structures have additional information like number and description
  return true;
}

/**
 * This implementation is simply transforming the char ** representation of the sample names into a 
 * contiguous vector<string>. As efficient as it gets.
 */
vector<string> VariantHeader::samples() const {
  return utils::hts_string_array_to_vector(m_header->samples, uint32_t(bcf_hdr_nsamples(m_header)));
}

vector<string> VariantHeader::chromosomes() const {
  return find_fields_of_type(m_header.get(), BCF_HL_CTG);
}

uint32_t VariantHeader::n_chromosomes() const {
  return count_fields_of_type(m_header.get(), BCF_HL_CTG);
}

vector<string> VariantHeader::filters() const {
  return find_fields_of_type(m_header.get(), BCF_HL_FLT);
}

uint32_t VariantHeader::n_filters() const {
  return count_fields_of_type(m_header.get(), BCF_HL_FLT);
}

vector<string> VariantHeader::shared_fields() const {
  return find_fields_of_type(m_header.get(), BCF_HL_INFO);
}

uint32_t VariantHeader::n_shared_fields() const {
  return count_fields_of_type(m_header.get(), BCF_HL_INFO);
}

vector<string> VariantHeader::individual_fields() const {
  return find_fields_of_type(m_header.get(), BCF_HL_FMT);
}

uint32_t VariantHeader::n_individual_fields() const {
  return count_fields_of_type(m_header.get(), BCF_HL_FMT);
}

vector<bcf_hrec_t*> VariantHeader::advanced_all_header_fields() const {
  return vector<bcf_hrec_t*>(m_header->hrec, m_header->hrec + m_header->nhrec);
}

}
