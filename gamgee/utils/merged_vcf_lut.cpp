#include "merged_vcf_lut.h"
#include "variant_utils.h"
#include "hts_memory.h"

using namespace std;

namespace gamgee
{

  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::MergedVCFLUTBase()
  {
    m_num_input_vcfs = 0u;
    m_num_merged_fields = 0u;
    clear();
  }

  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::MergedVCFLUTBase(unsigned numInputGVCFs, unsigned numMergedFields)
  {
    m_num_input_vcfs = numInputGVCFs;
    m_num_merged_fields = numMergedFields;
    clear();
    resize_inputs_2_merged_lut_if_needed(numInputGVCFs, numMergedFields);
    resize_merged_2_inputs_lut_if_needed(numInputGVCFs, numMergedFields);
  }

  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  void MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::clear()
  {
    for(auto& vec : m_inputs_2_merged_lut)
      vec.clear();
    m_inputs_2_merged_lut.clear();
    for(auto& vec : m_merged_2_inputs_lut)
      vec.clear();
    m_merged_2_inputs_lut.clear();
  }

  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  void MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::reset_vector(vector<int>& vec, unsigned from)
  {
    for(auto i=from;i<vec.size();++i)
      vec[i] = MergedVCFLUTBase::m_lut_missing_value;
  }

  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  void MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::resize_and_reset_vector(vector<int>& vec, unsigned new_size)
  {
    auto old_size = vec.size();
    if(new_size > old_size)
    {
      vec.resize(new_size);
      reset_vector(vec, old_size);
    }
  }

  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  void MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::resize_and_reset_lut
  (vector<vector<int>>& lut, unsigned new_lut_size, unsigned new_vector_size, unsigned& numRowsVar, unsigned& numColsVar)
  {
    auto old_lut_size = lut.size();
    if(new_lut_size > old_lut_size)
    {
      lut.resize(new_lut_size);
      numRowsVar = new_lut_size;
    }
    auto old_vector_size = (lut.size() > 0u) ? lut[0].size() : 0u;
    //Begin resizing of vectors at start_idx
    auto start_idx = old_lut_size;
    if(new_vector_size > old_vector_size)	//every vector needs to be resized
    {
      start_idx = 0u;
      numColsVar = new_vector_size;
    }
    else
      new_vector_size = old_vector_size;	//new vector size is smaller, don't bother reducing the size of existing rows
    for(auto i=start_idx;i<new_lut_size;++i)
      resize_and_reset_vector(lut[i], new_vector_size);
  }
  //explicit initialization to avoid link errors
  template class MergedVCFLUTBase<true,true>;

  //explicit initialization to avoid link errors
  template class MergedVCFAllelesIdxLUT<true,true>;

  //MergedVCFHeaderFieldsIdxLUT functions
  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  MergedVCFHeaderFieldsIdxLUT<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_header_fields_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx)
  {
    ASSERT(m_merged_vcf_header_ptr);
    for(auto j=0;j<curr_header->n[BCF_DT_ID];++j)
    {
      auto curr_id = &(curr_header->id[BCF_DT_ID][j]);
      //id has been deleted - ignore
      if(!bcf_hdr_idinfo_exists(curr_header, BCF_HL_FLT, j) && !bcf_hdr_idinfo_exists(curr_header, BCF_HL_INFO, j)
	  && !bcf_hdr_idinfo_exists(curr_header, BCF_HL_FMT, j))
	continue;
      int col_type = bcf_hdr_idinfo_exists(curr_header,BCF_HL_FLT, j) ? BCF_HL_FLT :
	bcf_hdr_idinfo_exists(curr_header,BCF_HL_INFO, j) ? BCF_HL_INFO :
	bcf_hdr_idinfo_exists(curr_header,BCF_HL_FMT, j) ? BCF_HL_FMT : -1;
      ASSERT(col_type != -1);
      bcf_hrec_t* hrec = bcf_hdr_id2hrec(curr_header, BCF_DT_ID, col_type, j);
      if(hrec) //not deleted
      {
	const char* key = curr_id->key;
	auto merged_idx = bcf_hdr_id2int(m_merged_vcf_header_ptr, BCF_DT_ID, key);
	ASSERT(merged_idx >= 0 && merged_idx < m_merged_vcf_header_ptr->n[BCF_DT_ID]);
	ASSERT(bcf_hdr_idinfo_exists(m_merged_vcf_header_ptr, col_type, merged_idx));
	m_header_fields_LUT.add_input_merged_idx_pair(input_vcf_idx, j, merged_idx);
      }
    }
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  MergedVCFHeaderFieldsIdxLUT<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_samples_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx)
  {
    for(auto j=0;j<bcf_hdr_nsamples(curr_header);++j)
    {
      if(m_sample2idx_merged.find(curr_header->samples[j]) == m_sample2idx_merged.end())
      {
	auto curr_size = m_sample2idx_merged.size();
	m_sample2idx_merged[curr_header->samples[j]] = curr_size;
      }
      m_samples_LUT.add_input_merged_idx_pair(input_vcf_idx, j, m_sample2idx_merged[curr_header->samples[j]]);
    }
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  MergedVCFHeaderFieldsIdxLUT<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  add_header(shared_ptr<bcf_hdr_t>& merged_hdr_ptr, const shared_ptr<bcf_hdr_t>& header_ptr)
  {
    if(merged_hdr_ptr)
      merge_variant_headers(merged_hdr_ptr, header_ptr);
    else
    {
      merged_hdr_ptr = utils::make_shared_variant_header(utils::variant_header_deep_copy(header_ptr.get()));
      m_merged_vcf_header_ptr = merged_hdr_ptr.get();
    }
    auto header_raw_ptr = header_ptr.get();
    ASSERT(m_merged_vcf_header_ptr);
    m_input_vcf_headers.push_back(header_raw_ptr);
    resize_luts_if_needed();

    unsigned input_vcf_idx = m_input_vcf_headers.size()-1;
    add_header_fields_mapping(header_raw_ptr, input_vcf_idx);
    add_samples_mapping(header_raw_ptr, input_vcf_idx);
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  MergedVCFHeaderFieldsIdxLUT<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  resize_luts_if_needed()
  {
    if(m_num_merged_fields_allocated < static_cast<unsigned>(m_merged_vcf_header_ptr->n[BCF_DT_ID]))
      m_num_merged_fields_allocated = m_merged_vcf_header_ptr->n[BCF_DT_ID] + 20;	//why 20, no particular reason
    if(m_input_vcf_headers.size() < m_num_input_vcfs_allocated)
      m_num_input_vcfs_allocated = 2*m_input_vcf_headers.size();
    if(m_num_merged_samples_allocated < static_cast<unsigned>(bcf_hdr_nsamples(m_merged_vcf_header_ptr)))
      m_num_merged_samples_allocated = 2*bcf_hdr_nsamples(m_merged_vcf_header_ptr);
    m_header_fields_LUT.resize_luts_if_needed(m_num_input_vcfs_allocated, m_num_merged_fields_allocated);
    m_samples_LUT.resize_luts_if_needed(m_num_input_vcfs_allocated, m_num_merged_samples_allocated);
  }

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  void 
  MergedVCFHeaderFieldsIdxLUT<fields_forward_LUT_ordering, fields_reverse_LUT_ordering, samples_forward_LUT_ordering, samples_reverse_LUT_ordering>::
  store_merged_field_idx_for_enum(const string& field_name, unsigned field_enum)
  {
    if(field_enum >= m_merged_field_idx_enum_vector.size())
      m_merged_field_idx_enum_vector.resize(2*field_enum + 1);
    auto val = -1;
    val = bcf_hdr_id2int(m_merged_vcf_header_ptr, BCF_HL_FLT, field_name.c_str());
    if(val == -1)
    {
      val = bcf_hdr_id2int(m_merged_vcf_header_ptr, BCF_HL_INFO, field_name.c_str());
      if(val == -1)
	val = bcf_hdr_id2int(m_merged_vcf_header_ptr, BCF_HL_FMT, field_name.c_str());
    }
    m_merged_field_idx_enum_vector[field_enum] = val; 
    if(m_merged_field_idx_enum_vector[field_enum] == -1)
      m_merged_field_idx_enum_vector[field_enum] = 
	MergedVCFLUTBase<fields_forward_LUT_ordering, fields_reverse_LUT_ordering>::m_lut_missing_value;
  }

  //explicit initialization to avoid link errors
  template class MergedVCFHeaderFieldsIdxLUT<true, true, true, true>;
}
