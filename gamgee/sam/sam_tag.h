#ifndef gamgee__sam_tag__guard
#define gamgee__sam_tag__guard

#include <string>

namespace gamgee {

/**
 * @brief class to represent a Sam TAG:TYPE:VALUE entry
 */
template<class TAG_TYPE>
class SamTag {
 public:
  explicit SamTag(const std::string& name, const TAG_TYPE& value, const bool missing = false) :
    m_name { name },
    m_value { value },
    m_missing { missing }
  {}

  explicit SamTag(const std::string& name, TAG_TYPE&& value, const bool missing = false) :
    m_name { name },
    m_value { std::move(value) },
    m_missing { missing }
  {}

  SamTag(const SamTag& other) = default;
  SamTag(SamTag&& other) = default;
  SamTag& operator=(const SamTag& other) = default;
  SamTag& operator=(SamTag&& other) = default;
  ~SamTag() = default;

  std::string name() const { return m_name; }
  TAG_TYPE value() const { return m_value; }
  bool missing() const { return m_missing; }

 private:
  std::string m_name;
  TAG_TYPE m_value;
  bool m_missing;
};

}

#endif // gamgee__sam_tag__guard
