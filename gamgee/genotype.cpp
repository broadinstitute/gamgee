#include "genotype.h"

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
  for (auto i = 0; i < count; i++) {
    if (allele_key(i) != other.allele_key(i)) {
      return false;
    }
  }
  return true;
}

bool Genotype::is_het() const {
  if (size() != 2) {
    return false;
  }
  const auto allele_1 = allele_key(0);
  const auto allele_2 = allele_key(1);
  return (allele_1 != allele_2);
}

bool Genotype::is_hom_var() const {
  if (size() != 2) {
    return false;
  }
  const auto allele_1 = allele_key(0);
  if (allele_1 == 0)
    return false;
  const auto allele_2 = allele_key(1);
  return allele_1 == allele_2;
}

bool Genotype::is_hom_ref() const {
  if (size() != 2) {
    return false;
  }
  const auto allele_1 = allele_key(0);
  if (allele_1 != 0)
    return false;
  const auto allele_2 = allele_key(1);
  return allele_2 == 0;
}

uint32_t Genotype::fast_diploid_key_generation() const {
  return (allele_key(0) << 16) | allele_key(1);
}

bool Genotype::is_missing() const {
  const auto count = size();
  for (auto i = 0; i < count; i++) {
    if (!utils::allele_missing(m_format_ptr, m_data_ptr, i)) {
      return false;
    }
  }
  return true;
}

std::vector<std::string> Genotype::alleles_strings() const {
  return utils::allele_strings(m_body, m_format_ptr, m_data_ptr);
}

std::vector<int32_t> Genotype::alleles_keys() const {
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

}
