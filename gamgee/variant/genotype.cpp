#include "genotype.h"

#include <algorithm>

#include "boost/algorithm/cxx11/all_of.hpp"

namespace gamgee {

using namespace std;

/**
 * @note Most of the patterns here are directly copied or adapted from htslib.
 */
Genotype::Genotype(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr) :
  m_body {body},
  m_format_ptr {format_ptr},
  m_data_ptr {data_ptr} {
}

bool Genotype::operator!=(const Genotype& other) const {
  return !(*this == other);
}

bool Genotype::operator==(const Genotype& other) const {
  const auto count = size();
  if (count != other.size())
    return false;
  for (auto i = 0u; i != count; ++i) {
    if (allele_key(i) != other.allele_key(i)) {
      return false;
    }
  }
  return true;
}

bool Genotype::het() const {
  if (size() != 2) {
    return false;
  }
  const auto allele_1 = allele_key(0);
  const auto allele_2 = allele_key(1);
  return (allele_1 != allele_2);
}

bool Genotype::non_ref_het() const {
  if (size() != 2) {
    return false;
  }
  const auto allele_1 = allele_key(0);
  if (allele_1 == 0)
    return false;
  const auto allele_2 = allele_key(1);
  if (allele_2 == 0)
    return false;
  return (allele_1 != allele_2);
}

bool Genotype::hom_var() const {
  const auto keys = allele_keys();
  const auto allele_1 = keys[0];
  return allele_1 != 0 && boost::algorithm::all_of_equal(keys, allele_1);
}

bool Genotype::hom_ref() const {
  const auto keys = allele_keys();
  return boost::algorithm::all_of_equal(keys, 0);
}

uint32_t Genotype::fast_diploid_key_generation() const {
  return (allele_key(0) << 16) | allele_key(1);
}

bool Genotype::missing() const {
  const auto count = size();
  for (auto i = 0u; i != count; ++i) {
    if (!utils::allele_missing(m_format_ptr, m_data_ptr, i)) {
      return false;
    }
  }
  return true;
}

std::vector<std::string> Genotype::allele_strings() const {
  return utils::allele_strings(m_body, m_format_ptr, m_data_ptr);
}

std::vector<int32_t> Genotype::allele_keys() const {
  return utils::allele_keys(m_body, m_format_ptr, m_data_ptr);
}

std::string Genotype::allele_string(const uint32_t index) const {
  return utils::allele_key_to_string(m_body, allele_key(index));
}

int32_t Genotype::allele_key(const uint32_t index) const {
  return utils::allele_key(m_format_ptr, m_data_ptr, index);
}

int32_t Genotype::operator[](const uint32_t index) const {
  return allele_key(index);
}

uint32_t Genotype::size() const {
  return utils::allele_count(m_format_ptr);
}

bool Genotype::snp(const AlleleMask& mask) const {
  const auto keys = allele_keys();
  return allele_is_type_or_ref(AlleleType::SNP, keys, mask);
}

bool Genotype::insertion(const AlleleMask& mask) const {
  const auto keys = allele_keys();
  return allele_is_type_or_ref(AlleleType::INSERTION, keys, mask);
}

bool Genotype::deletion(const AlleleMask& mask) const {
  const auto keys = allele_keys();
  return allele_is_type_or_ref(AlleleType::DELETION, keys, mask);
}

bool Genotype::indel(const AlleleMask& mask) const {
  const auto keys = allele_keys();
  auto found_type = false;
  const auto result = std::all_of(keys.begin(), keys.end(), [&mask,&found_type](const int32_t& k) {
      if (mask[k] == AlleleType::INSERTION || mask[k] == AlleleType::DELETION) {
        found_type = true;
        return true;
      }
      return mask[k] == AlleleType::REFERENCE ? true : false;
  });
  return result && found_type;
}

bool Genotype::biallelic() const {
  const auto keys = allele_keys();
  const auto first_it = find_if(keys.begin(), keys.end(), [](const int32_t& k){return k != 0;}); // 0 = reference key
  if (first_it == keys.end()) 
    return true;
  return all_of(first_it, keys.end(), [&first_it](const int32_t& k){return k == 0 || k == *first_it;});
}

bool Genotype::mixed() const {
  const auto keys = allele_keys();
  const auto first_it = find_if(keys.begin(), keys.end(), [](const int32_t& k){return k != 0;}); // 0 = reference key
  if (first_it == keys.end()) 
    return false;
  return keys.end() != find_if(first_it, keys.end(), [&first_it](const int32_t& k){return k != 0 && k != *first_it;});
}

bool Genotype::allele_is_type_or_ref(const AlleleType& type, const std::vector<int32_t>& keys, const AlleleMask& mask) const {
  auto found_type = false;
  const auto result = std::all_of(keys.begin(), keys.end(), [&type,&mask,&found_type](const int32_t& k) {
      if (mask[k] == type) {
        found_type = true;
        return true;
      }
      return mask[k] == AlleleType::REFERENCE ? true : false;
  });
  return result && found_type;
}

} // end of namespace
