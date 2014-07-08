#ifndef __gamgee__genotype_field_value__
#define __gamgee__genotype_field_value__

#include "genotype_field_value_iterator.h"
#include "diploid_pl_genotype.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"
#include "utils/format_field_type.h"

#include "htslib/vcf.h"

#include <memory>

namespace gamgee {

class GenotypeFieldValue {

  public:
    GenotypeFieldValue(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, uint8_t* const data_ptr) :
      m_body {body},
      m_format_ptr {const_cast<bcf_fmt_t*>(format_ptr)},
      m_data_ptr {data_ptr},
      m_num_bytes {utils::size_for_type(static_cast<utils::FormatFieldType>(format_ptr->type))}
    {
    }

    GenotypeFieldValue(const GenotypeFieldValue& other) = delete;

    GenotypeFieldValue(GenotypeFieldValue&& other) :
      m_body {std::move(other.m_body)},
      m_format_ptr {other.m_format_ptr},
      m_data_ptr {other.m_data_ptr},
      m_num_bytes {other.m_num_bytes}
    {
    }

    GenotypeFieldValue& operator=(const GenotypeFieldValue& other) = delete;

    GenotypeFieldValue& operator=(GenotypeFieldValue&& other) {
      if (this != &other)
        return *this;
      m_body = std::move(other.m_body);
      m_format_ptr = other.m_format_ptr;
      m_data_ptr = other.m_data_ptr;
      m_num_bytes = other.m_num_bytes;
      return *this;
    }

    int32_t operator[](const uint32_t index) const {
      utils::check_boundaries(index, m_format_ptr->n - 1);
      return convert_from_byte_array(index);
    }

    GenotypeFieldValueIterator begin() const {
      return GenotypeFieldValueIterator{m_body, m_data_ptr, m_num_bytes, static_cast<utils::FormatFieldType>(m_format_ptr->type)};
    }

    GenotypeFieldValueIterator end() const {
      return GenotypeFieldValueIterator{m_body, m_data_ptr + m_format_ptr->size, m_num_bytes, static_cast<utils::FormatFieldType>(m_format_ptr->type)};
    }

    /**
     * @brief whether or not the sample is hom_ref
     */
    inline bool is_hom_ref() const {
      return is_this_genotype(DiploidPLGenotype::HOM_REF);
    }

    /**
     * @brief whether or not the sample is het
     */
    inline bool is_het() const {
      return is_this_genotype(DiploidPLGenotype::HET);
    }

    /**
     * @brief whether or not the sample is hom_var
     */
    inline bool is_hom_var() const {
      return is_this_genotype(DiploidPLGenotype::HOM_VAR);
    }

  private:
    std::shared_ptr<bcf1_t> m_body;
    bcf_fmt_t* m_format_ptr;
    uint8_t* m_data_ptr;
    uint8_t m_num_bytes;

    bool is_this_genotype(const DiploidPLGenotype& genotype) const {
      return (*this)[static_cast<int32_t>(genotype)] == 0;
    }

    inline int32_t convert_from_byte_array(int index) const {
      return utils::convert_data_to_integer(m_data_ptr, index, m_num_bytes, static_cast<utils::FormatFieldType>(m_format_ptr->type));
    }
};

}

#endif
