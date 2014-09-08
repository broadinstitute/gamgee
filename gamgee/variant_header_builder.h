#ifndef gamgee__variant_header_builder__guard
#define gamgee__variant_header_builder__guard

#include "variant_header.h"

#include "utils/utils.h"
#include "utils/hts_memory.h"

#include "htslib/vcf.h"

#include <memory>

namespace gamgee {

/**
 * @brief Utility class to build VariantHeader objects from scratch
 */
class VariantHeaderBuilder {
 public: 

  /**
   * @brief         initializes a variant header builder 
   * @note          VariantReader and VariantWriter will take care of creating objects for reading and writing respectively if your goal is to write a VCF/BCF file to disk or stream.
   * @warning       you should only call this explicitly if you are building a VariantHeader from scratch.
   */
  VariantHeaderBuilder() noexcept;

  VariantHeaderBuilder(const VariantHeaderBuilder& other) = delete;
  VariantHeaderBuilder& operator=(const VariantHeaderBuilder& other) = delete;

  VariantHeaderBuilder(VariantHeaderBuilder&& other) = default;
  VariantHeaderBuilder& operator=(VariantHeaderBuilder&& other) = default;

  VariantHeaderBuilder& add_chromosome(const std::string& id, const std::string& length = "", const std::string& url = "", const std::string& extra = "");
  VariantHeaderBuilder& add_filter(const std::string& id, const std::string& description = "", const std::string& extra = "");
  VariantHeaderBuilder& add_shared_field(const std::string& id, const std::string& number, const std::string& type, const std::string& description = "", const std::string& source = "", const std::string& version = "", const std::string& extra = "");
  VariantHeaderBuilder& add_individual_field(const std::string& id, const std::string& number, const std::string& type, const std::string& description = "", const std::string& extra = "");
  VariantHeaderBuilder& add_sample(const std::string& sample);
  VariantHeaderBuilder& add_source(const std::string& source);

  VariantHeaderBuilder& advanced_add_arbitrary_line(const std::string& line);

  /**
   * @brief create a new VariantHeader by copying the contents of this builder into a new object.  Allows for reuse.
   */
  VariantHeader build() const {
    return VariantHeader{utils::make_shared_variant_header(utils::variant_header_deep_copy(m_header.get()))};
  }

  /**
   * @brief create a new VariantHeader by moving the contents of this builder into a new object.
   * More efficient but does not allow for reuse.
   */
  VariantHeader one_time_build() const { return VariantHeader{move(m_header)}; }

 private:
  std::shared_ptr<bcf_hdr_t> m_header;

};

}

#endif // gamgee__variant_header_builder__guard
