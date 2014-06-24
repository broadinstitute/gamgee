#ifndef __gamgee__sam_tag__
#define __gamgee__sam_tag__

#include <string>

namespace gamgee {

/**
 * @brief class to represent a Sam TAG:TYPE:VALUE entry
 */
template<class TAG_TYPE>
class SamTag {
 public:
  explicit SamTag(const std::string& name, const TAG_TYPE& value, const bool is_present = true) :
    m_name { name },
    m_value { value },
    m_is_present { is_present }
  {}

  explicit SamTag(std::string& name, TAG_TYPE&& value, const bool is_present = true) :
    m_name { name },
    m_value { std::move(value) },
    m_is_present { is_present }
  {}

  SamTag(const SamTag& other) = default;
  SamTag& operator=(const SamTag& other) = default;
  ~SamTag() = default;

  SamTag(SamTag&& other) noexcept :
    m_name { std::move(other.m_name) },
    m_value { std::move(other.m_value) },
    m_is_present { other.m_is_present }
  {}

  SamTag& operator=(SamTag&& other) noexcept {
    if ( &other == this )
      return *this;

    m_name = std::move(other.m_name);
    m_value = std::move(other.m_value);
    m_is_present = other.m_is_present;
    return *this;
  }

  std::string name() const { return m_name; }
  TAG_TYPE value() const { return m_value; }
  bool is_present() const { return m_is_present; }

 private:
  std::string m_name;
  TAG_TYPE m_value;
  bool m_is_present;
};

}

#endif // __gamgee__sam_tag__
