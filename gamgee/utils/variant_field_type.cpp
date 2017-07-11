#include "variant_field_type.h"

#include "htslib/vcf.h"

#include <string>
#include <stdexcept>

namespace gamgee {
namespace utils {

int32_t convert_data_to_integer(const uint8_t* data_ptr, const int index, const uint32_t num_bytes_per_value, const VariantFieldType& type) {
  auto return_val = -1;
  const auto p = data_ptr + (index * num_bytes_per_value);
  auto end_val = -1;
  switch (type) {
    case VariantFieldType::INT8:
      return_val = *(reinterpret_cast<const int8_t*>(p));
      end_val = bcf_int8_vector_end;
      break;
    case VariantFieldType::INT16:
      return_val = *(reinterpret_cast<const int16_t*>(p));
      end_val = bcf_int16_vector_end;
      break;
    case VariantFieldType::INT32:
      return *(reinterpret_cast<const int32_t*>(p));
    case VariantFieldType::FLOAT:
      return *(reinterpret_cast<const float*>(p));
    case VariantFieldType::STRING:
      throw std::invalid_argument("user requested an integer value but underlying type is actually a string");
    default:
      return 0; // undefined or NULL type -- impossible to happen just so the compiler doesn't warn us 
  }
  //diff = 0 if return_val == vector_end, 1 if return_val == missing, < 0 if return_val is 'normal'
  auto diff = end_val - return_val;
  if(diff >= 0)
    return bcf_int32_vector_end - diff;
  else
    return return_val;
}

float convert_data_to_float(const uint8_t* data_ptr, const int index, const uint32_t num_bytes_per_value, const VariantFieldType& type) {
  auto return_val = -1;
  auto end_val = -1;
  const auto p = data_ptr + (index * num_bytes_per_value);
  switch (type) {
    case VariantFieldType::INT8:
      return_val = *(reinterpret_cast<const int8_t*>(p)); 
      end_val = bcf_int8_vector_end;
      break;
    case VariantFieldType::INT16:
      return_val = *(reinterpret_cast<const int16_t*>(p));
      end_val = bcf_int16_vector_end;
      break;
    case VariantFieldType::INT32:
      return_val = *(reinterpret_cast<const int32_t*>(p)); 
      end_val = bcf_int32_vector_end;
      break;
    case VariantFieldType::FLOAT:
      return *(reinterpret_cast<const float*>(p)); 
    case VariantFieldType::STRING:
      throw std::invalid_argument("user requested a float value but underlying type is actually a string");
    default:
      return 0.0f; // undefined or NULL type -- impossible to happen just so the compiler doesn't warn us 
  }
  if(end_val >= return_val)	//can't use diff because INT32 types are also checked here and diff will wrap-around 
  {
    auto tmp = 0.0f;	//pack bcf_float_missing or bcf_float_vector_end into tmp and return tmp
    //diff = 0 if return_val == vector_end, 1 if return_val == missing, 
    auto diff = end_val - return_val;
    bcf_float_set(&tmp, bcf_float_vector_end - diff);
    return tmp;
  }
  else
    return return_val;
}

std::string convert_data_to_string(const uint8_t* data_ptr, const int index, const uint32_t num_bytes_per_value, const VariantFieldType& type) {
  auto result = std::string{};
  const auto p = data_ptr + (index * num_bytes_per_value);
  switch (type) {
    case VariantFieldType::INT8:
      return std::to_string(*(reinterpret_cast<const int8_t*>(p))); // convert from integer to string if necessary
    case VariantFieldType::INT16:
      return std::to_string(*(reinterpret_cast<const int16_t*>(p))); // convert from integer to string if necessary
    case VariantFieldType::INT32:
      return std::to_string(*(reinterpret_cast<const int32_t*>(p))); // convert from integer to string if necessary
    case VariantFieldType::FLOAT:
      return std::to_string(*(reinterpret_cast<const float*>(p))); // convert from float to string if necessary
    case VariantFieldType::STRING:
      for (auto i = 0u; i != num_bytes_per_value; ++i) {
        const auto c = reinterpret_cast<const char *>(p)[i];
        if (c == bcf_str_vector_end) 
          break;
        result += c;
      }
      return result;
    default:
      return std::string{}; // undefined or NULL type -- impossible to happen just so the compiler doesn't warn us 
  }
}

uint32_t size_for_type(const VariantFieldType& type, const bcf_fmt_t* const format_ptr) {
  switch (type) {
    case VariantFieldType::NIL:
    case VariantFieldType::INT8: 
    case VariantFieldType::INT16:
      return static_cast<uint8_t>(type);
    case VariantFieldType::INT32:
    case VariantFieldType::FLOAT:
      return 4;
    case VariantFieldType::STRING:
      return format_ptr->n; // htslib keeps the number of bytes in the string in the n member variable. This is weird, but that's how it is.
    default: 
      return 0; 
  }
}

uint32_t size_for_type(const VariantFieldType& type, const bcf_info_t* const info_ptr) {
  switch (type) {
    case VariantFieldType::NIL:
    case VariantFieldType::INT8: 
    case VariantFieldType::INT16:
      return static_cast<uint8_t>(type);
    case VariantFieldType::INT32:
    case VariantFieldType::FLOAT:
      return 4;
    case VariantFieldType::STRING:
      return info_ptr->len; // htslib keeps the number of bytes in the string in the n member variable. This is weird, but that's how it is.
    default: 
      return 0; 
  }
}

} // end namespace utils
} // end namespace gamgee

