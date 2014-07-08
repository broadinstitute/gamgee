#ifndef __gamgee__genotype_field_value__
#define __gamgee__genotype_field_value__

#include "diploid_pl_genotype.h"
#include "utils/hts_memory.h"
#include "utils/utils.h"
#include "utils/format_field_type.h"

#include "htslib/vcf.h"

#include <memory>

namespace gamgee {

class GenotypeFieldValue : public VariantFieldValue<int32_t> {

  public:
    GenotypeFieldValue(const std::shared_ptr<bcf1_t>& body, const bcf_fmt_t* const format_ptr, uint8_t* const data_ptr) :
      VariantFieldValue<int32_t>(body, format_ptr, data_ptr)
    {
    }

    GenotypeFieldValue(const GenotypeFieldValue& other) = delete;

    GenotypeFieldValue(GenotypeFieldValue&& other) :
      VariantFieldValue<int32_t>(std::move(other))
    {
    }

    GenotypeFieldValue& operator=(const GenotypeFieldValue& other) = delete;

    GenotypeFieldValue& operator=(GenotypeFieldValue&& other) {
      if (this != &other)
        return *this;
      VariantFieldValue<int32_t>::operator=(std::move(other));
      return *this;
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

    bool is_this_genotype(const DiploidPLGenotype& genotype) const {
      return (*this)[static_cast<int32_t>(genotype)] == 0;
    }
};

}

#endif
