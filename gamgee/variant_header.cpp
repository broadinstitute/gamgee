#include "variant_header.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"
#include "htslib/vcf.h"
#include "missing.h"

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

vector<string> VariantHeader::chromosomes() const {
  return find_fields_of_type(m_header.get(), BCF_HL_CTG);
}

vector<string> VariantHeader::filters() const {
  return find_fields_of_type(m_header.get(), BCF_HL_FLT);
}

vector<string> VariantHeader::shared_fields() const {
  return find_fields_of_type(m_header.get(), BCF_HL_INFO);
}

vector<string> VariantHeader::individual_fields() const {
  return find_fields_of_type(m_header.get(), BCF_HL_FMT);
}

bool VariantHeader::has_individual_field(const string field) const {
  const auto& fields = individual_fields();
  return find(fields.begin(), fields.end(), field) != fields.end();
}

int32_t VariantHeader::individual_field_index(const std::string key) const {
  const auto hdr = m_header.get();
  const auto id = bcf_hdr_id2int(hdr, BCF_DT_ID, key.c_str());
  if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,id) ) return missing_values::int32;   // no such FMT field in the header
  return id;
}

int32_t VariantHeader::shared_field_index(const std::string tag) const {
  const auto hdr = m_header.get();
  const auto tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag.c_str());
  if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,tag_id) ) return missing_values::int32;    // no such INFO field in the header
  return tag_id;
}

}
