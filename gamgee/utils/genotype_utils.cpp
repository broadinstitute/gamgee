#include "genotype_utils.h"

namespace gamgee {

namespace utils {

using namespace std;

// Hidden template declarations, used only in this cpp file.

template<class TYPE>
bool allele_missing(const uint8_t* data_ptr, const uint32_t allele_index, const TYPE missing);

template<class TYPE>
vector<int32_t> allele_keys(const std::shared_ptr<bcf1_t>& body,
    const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr,
    const TYPE missing, const TYPE vector_end);

template<class TYPE>
vector<string> allele_strings(const std::shared_ptr<bcf1_t>& body,
    const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr,
    const TYPE missing, const TYPE vector_end);

bool allele_missing(const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr, const uint32_t allele_index) {
  switch (format_ptr->type) {
  case BCF_BT_INT8:
    return allele_missing<int8_t>(data_ptr, allele_index, bcf_int8_missing);
  case BCF_BT_INT16:
    return allele_missing<int16_t>(data_ptr, allele_index, bcf_int16_missing);
  case BCF_BT_INT32:
    return allele_missing<int32_t>(data_ptr, allele_index, bcf_int32_missing);
  default:
    throw invalid_argument("unknown GT field type: " + to_string(format_ptr->type));
  }
}

template<class TYPE>
bool allele_missing(const uint8_t* data_ptr, const uint32_t allele_index, const TYPE missing) {
  // mostly copied from htslib
  const auto p = reinterpret_cast<const TYPE*>(data_ptr);
  return !(p[allele_index]>>1) || p[allele_index] == missing;
}

vector<int32_t> allele_keys(const std::shared_ptr<bcf1_t>& body,
    const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr) {
  switch (format_ptr->type) {
  case BCF_BT_INT8:
    return allele_keys<int8_t>(body, format_ptr, data_ptr, bcf_int8_missing, bcf_int8_vector_end);
  case BCF_BT_INT16:
    return allele_keys<int16_t>(body, format_ptr, data_ptr, bcf_int16_missing, bcf_int16_vector_end);
  case BCF_BT_INT32:
    return allele_keys<int32_t>(body, format_ptr, data_ptr, bcf_int32_missing, bcf_int32_vector_end);
  default:
    throw invalid_argument("unknown GT field type: " + to_string(format_ptr->type));
  }
}

template<class TYPE>
vector<int32_t> allele_keys(const std::shared_ptr<bcf1_t>& body,
    const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr,
    const TYPE missing, const TYPE vector_end) {
  const auto p = reinterpret_cast<const TYPE*>(data_ptr);
  // mostly copied from htslib
  auto results = vector<int32_t>{};
  results.reserve(format_ptr->n);
  for (int ial=0; ial<format_ptr->size; ial++)
  {
      if ( p[ial]==vector_end ) break; /* smaller ploidy */
      if ( !(p[ial]>>1) || p[ial]==missing ) {
        /* missing allele */
        results.push_back(missing_values::int32);
        continue;
      }
      if ( (p[ial]>>1)-1 >= body->n_allele ) {
        results.push_back(missing_values::int32);
        continue;
      }
      results.push_back((p[ial]>>1)-1);
  }
  results.shrink_to_fit();
  return results;
}

vector<string> allele_strings(const std::shared_ptr<bcf1_t>& body,
    const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr) {
  switch (format_ptr->type) {
  case BCF_BT_INT8:
    return allele_strings<int8_t>(body, format_ptr, data_ptr, bcf_int8_missing, bcf_int8_vector_end);
  case BCF_BT_INT16:
    return allele_strings<int16_t>(body, format_ptr, data_ptr, bcf_int16_missing, bcf_int16_vector_end);
  case BCF_BT_INT32:
    return allele_strings<int32_t>(body, format_ptr, data_ptr, bcf_int32_missing, bcf_int32_vector_end);
  default:
    throw invalid_argument("unknown GT field type: " + to_string(format_ptr->type));
  }
}

template<class TYPE>
vector<string> allele_strings(const std::shared_ptr<bcf1_t>& body,
    const bcf_fmt_t* const format_ptr, const uint8_t* data_ptr,
    const TYPE missing, const TYPE vector_end) {
  const auto p = reinterpret_cast<const TYPE*>(data_ptr);
  // mostly copied from htslib
  auto results = vector<string>{};
  results.reserve(format_ptr->n);
  bcf_unpack(body.get(), BCF_UN_STR);
  for (int ial=0; ial<format_ptr->size; ial++)
  {
      if ( p[ial]==vector_end ) break; /* smaller ploidy */
      if ( !(p[ial]>>1) || p[ial]==missing ) {
        /* missing allele */
        results.push_back(missing_values::string_empty);
        continue;
      }
      if ( (p[ial]>>1)-1 >= body->n_allele ) {
        results.push_back(missing_values::string_empty);
        continue;
      }
      results.push_back(allele_key_to_string(body, (p[ial]>>1)-1));
  }
  results.shrink_to_fit();
  return results;
}

string allele_key_to_string(
    const std::shared_ptr<bcf1_t>& body, const int32_t allele_int) {
  if (allele_int == missing_values::int32) {
    return missing_values::string_empty;
  }
  bcf_unpack(body.get(), BCF_UN_STR);
  return string{body->d.allele[allele_int]};
}

}

}
