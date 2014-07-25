#ifndef gamgee__variant_header_builder__guard
#define gamgee__variant_header_builder__guard

#include "utils/utils.h"
#include "variant_header.h"

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
  VariantHeaderBuilder(VariantHeaderBuilder&& other) noexcept;
  VariantHeaderBuilder(const VariantHeaderBuilder& other) = delete; ///< @brief explicitly forbid copying of a builder

  void add_chromosome(const std::string& id, const std::string& length = "", const std::string& url = "", const std::string& extra = "");
  void add_filter(const std::string& id, const std::string& description = "", const std::string& extra = "");
  void add_info_field(const std::string& id, const std::string& number, const std::string& type, const std::string& description = "", const std::string& source = "", const std::string& version = "", const std::string& extra = "");
  void add_format_field(const std::string& id, const std::string& number, const std::string& type, const std::string& description = "", const std::string& extra = "");
  void add_sample(const std::string& sample);
  void add_source(const std::string& source);

  void advanced_add_arbitrary_line(const std::string& line);

  VariantHeader build() const { return VariantHeader{m_header}; }

 private:
  std::shared_ptr<bcf_hdr_t> m_header;

};

}

#endif // gamgee__variant_header_builder__guard
