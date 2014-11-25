#include "variant_header_merger.h"
#include "utils/variant_utils.h"
#include "utils/hts_memory.h"

using namespace std;

namespace gamgee
{
  //VariantHeaderMerger functions
  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  VariantHeaderMerger<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_header_fields_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx)
  {
    assert(m_merged_vcf_header_ptr);
    for(auto j=0;j<curr_header->n[BCF_DT_ID];++j)
    {
      auto curr_id = &(curr_header->id[BCF_DT_ID][j]);
      for(auto bcf_hl_type : { BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT })
      {
        //id has been deleted - ignore
        if(!bcf_hdr_idinfo_exists(curr_header, bcf_hl_type, j))
          continue;
        bcf_hrec_t* hrec = bcf_hdr_id2hrec(curr_header, BCF_DT_ID, bcf_hl_type, j);
        if(hrec) //not deleted
        {
          const char* key = curr_id->key;
          auto merged_idx = bcf_hdr_id2int(m_merged_vcf_header_ptr.get(), BCF_DT_ID, key);
          assert(merged_idx >= 0 && merged_idx < m_merged_vcf_header_ptr->n[BCF_DT_ID]);
          assert(bcf_hdr_idinfo_exists(m_merged_vcf_header_ptr, bcf_hl_type, merged_idx));
          m_header_fields_LUT.add_input_merged_idx_pair(input_vcf_idx, j, merged_idx);
        }
      }
    }
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  VariantHeaderMerger<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_samples_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx)
  {
    for(auto j=0;j<bcf_hdr_nsamples(curr_header);++j)
    {
      if(curr_header->samples[j] && bcf_hdr_id2int(curr_header, BCF_DT_SAMPLE, curr_header->samples[j]) >= 0)
      {
        if(m_sample2idx_merged.find(curr_header->samples[j]) == m_sample2idx_merged.end())
        {
          auto curr_size = m_sample2idx_merged.size();
          m_sample2idx_merged[curr_header->samples[j]] = curr_size;
        }
        m_samples_LUT.add_input_merged_idx_pair(input_vcf_idx, j, m_sample2idx_merged[curr_header->samples[j]]);
      }
      else
        m_samples_LUT.reset_merged_idx_for_input(input_vcf_idx, j);
    }
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  VariantHeaderMerger<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_header(const shared_ptr<bcf_hdr_t>& header_ptr)
  {
    if(m_merged_vcf_header_ptr)
      merge_variant_headers(m_merged_vcf_header_ptr, header_ptr);
    else
      m_merged_vcf_header_ptr = utils::make_shared_variant_header(utils::variant_header_deep_copy(header_ptr.get()));
    auto header_raw_ptr = header_ptr.get();
    assert(m_merged_vcf_header_ptr);
    m_input_vcf_headers.push_back(header_ptr);
    resize_luts_if_needed();

    unsigned input_vcf_idx = m_input_vcf_headers.size()-1;
    add_header_fields_mapping(header_raw_ptr, input_vcf_idx);
    add_samples_mapping(header_raw_ptr, input_vcf_idx);
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  VariantHeaderMerger<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_header(const VariantHeader& header)
  {
    add_header(header.m_header);
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  VariantHeaderMerger<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_headers(const vector<shared_ptr<bcf_hdr_t>>& headers)
  {
    for(const auto& header : headers)
      add_header(header);
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  VariantHeaderMerger<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_headers(const vector<VariantHeader>& headers)
  {
    for(const auto& header : headers)
      add_header(header);
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  VariantHeaderMerger<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  resize_luts_if_needed()
  {
    if(m_num_merged_fields_allocated < static_cast<unsigned>(m_merged_vcf_header_ptr->n[BCF_DT_ID]))
      m_num_merged_fields_allocated = m_merged_vcf_header_ptr->n[BCF_DT_ID] + 20;       //why 20, no particular reason
    if(m_num_input_vcfs_allocated < m_input_vcf_headers.size())
      m_num_input_vcfs_allocated = 2*m_input_vcf_headers.size();
    if(m_num_merged_samples_allocated < static_cast<unsigned>(bcf_hdr_nsamples(m_merged_vcf_header_ptr)))
      m_num_merged_samples_allocated = 2*bcf_hdr_nsamples(m_merged_vcf_header_ptr);
    m_header_fields_LUT.resize_luts_if_needed(m_num_input_vcfs_allocated, m_num_merged_fields_allocated);
    m_samples_LUT.resize_luts_if_needed(m_num_input_vcfs_allocated, m_num_merged_samples_allocated);
    m_merged_field_idx_enum_lut.resize_luts_if_needed(1u, m_num_merged_fields_allocated);
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  VariantHeaderMerger<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  store_merged_field_idx_for_enum(const string& field_name, unsigned field_enum)
  {
    if(field_enum >= m_num_enums_allocated)
    {
      m_num_enums_allocated = 2u*field_enum + 1u;
      m_merged_field_idx_enum_lut.resize_luts_if_needed(1u, m_num_enums_allocated);
    }
    auto val =  bcf_hdr_id2int(m_merged_vcf_header_ptr.get(), BCF_DT_ID, field_name.c_str());
    if(val == -1)
      m_merged_field_idx_enum_lut.reset_merged_idx_for_input(0u, field_enum);
    else
      m_merged_field_idx_enum_lut.add_input_merged_idx_pair(0u, field_enum, val);
  }

  //explicit initialization to avoid link errors
  template class VariantHeaderMerger<true, true, true, true>;
  template class VariantHeaderMerger<false, false, false, false>;
}
