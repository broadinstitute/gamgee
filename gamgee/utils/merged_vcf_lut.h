#ifndef __gamgee_merged_vcf_lut__
#define __gamgee_merged_vcf_lut__

#include "htslib/vcf.h"

#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>
#include <unordered_map>

#define KG_DEBUG 1
#ifdef KG_DEBUG
#include<assert.h>
#define ASSERT(x) assert(x)
#else
#define ASSERT(x)   ;
#endif

using namespace std;
namespace gamgee
{
  //friend class of MergedVCFLUTBase
  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  class MergedVCFHeaderFieldsIdxLUT;

  //LUT = Look Up Table (to avoid confusion with map, unordered_map etc)
  //Base class to store look up information between fields of merged header and input headers
  //Each MergedVCFLUTBase object contains 2 tables (vector of vector): one for mapping input field idx to merged field idx (m_inputs_2_merged_lut)
  //and the second for mapping merged field idx to input field idx (m_merged_2_inputs_lut)
  //The boolean template parameters specify how the 2 tables are laid out in memory - whether the outer vector corresponds to fields or input vcfs
  //For example, in object of type MergedVCFLUTBase<true, true>, both LUTs are laid out such that
  //m_inputs_2_merged_lut[0] contains mappings for all fields for input VCF file 0. This would lead to fast traversal of all fields for
  //a given input VCF (cache locality). However, traversing over all input VCFs for a given field would be slow (many cache misses).
  //The value of the template parameters depends on the application using the LUT.
  //
  //Almost all the 'complexity' of the code comes from being able to handle the different layouts in a transparent manner
  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  class MergedVCFLUTBase
  {
    template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
    friend class MergedVCFHeaderFieldsIdxLUT;
    public:
    inline void reset_luts()
    {
      for(auto& vec : m_inputs_2_merged_lut)
	reset_vector(vec);
      for(auto& vec : m_merged_2_inputs_lut)
	reset_vector(vec);
    }

    inline void add_input_merged_idx_pair(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    {
      set_merged_idx_for_input(inputGVCFIdx, inputIdx, mergedIdx);
      set_input_idx_for_merged(inputGVCFIdx, inputIdx, mergedIdx);
    }

    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    inline int get_input_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx)
    { return get_lut_value(m_merged_2_inputs_lut, inputGVCFIdx, mergedIdx); }
    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    inline int get_input_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx)
    { return get_lut_value(m_merged_2_inputs_lut, mergedIdx, inputGVCFIdx); }

    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    inline int get_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx)
    { return get_lut_value(m_inputs_2_merged_lut, inputGVCFIdx, inputIdx); }
    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    inline int get_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx)
    { return get_lut_value(m_inputs_2_merged_lut, inputIdx, inputGVCFIdx); }

    inline void reset_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx)
    {	set_merged_idx_for_input(inputGVCFIdx, inputIdx, MergedVCFLUTBase::m_lut_missing_value); }

    inline void reset_input_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx)
    {	set_input_idx_for_merged(inputGVCFIdx, MergedVCFLUTBase::m_lut_missing_value, mergedIdx); }

    static inline bool is_missing(const int lut_val) {  return (lut_val == MergedVCFLUTBase::m_lut_missing_value); }
    protected:
    static const int m_lut_missing_value = -1;
    MergedVCFLUTBase(); 
    //Only inherited classes should call constructor,destructor etc
    MergedVCFLUTBase(unsigned numInputGVCFs, unsigned numMergedFields);
    ~MergedVCFLUTBase() { clear(); }
    void clear();

    unsigned m_num_input_vcfs;
    unsigned m_num_merged_fields;

    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    void resize_inputs_2_merged_lut_if_needed(unsigned numInputGVCFs, unsigned numMergedFields)
    {	resize_and_reset_lut(m_inputs_2_merged_lut, numInputGVCFs, numMergedFields, m_num_input_vcfs, m_num_merged_fields); }

    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    void resize_inputs_2_merged_lut_if_needed(unsigned numInputGVCFs, unsigned numMergedFields)
    {	resize_and_reset_lut(m_inputs_2_merged_lut, numMergedFields, numInputGVCFs, m_num_merged_fields, m_num_input_vcfs); }

    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    void resize_merged_2_inputs_lut_if_needed(unsigned numInputGVCFs, unsigned numMergedFields)
    {	resize_and_reset_lut(m_merged_2_inputs_lut, numInputGVCFs, numMergedFields, m_num_input_vcfs, m_num_merged_fields); }

    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    void resize_merged_2_inputs_lut_if_needed(unsigned numInputGVCFs, unsigned numMergedFields)
    {	resize_and_reset_lut(m_merged_2_inputs_lut, numMergedFields, numInputGVCFs, m_num_merged_fields, m_num_input_vcfs); }

    void resize_luts_if_needed(unsigned numInputGVCFs, unsigned numMergedFields)
    {
      resize_merged_2_inputs_lut_if_needed(numInputGVCFs, numMergedFields);
      resize_inputs_2_merged_lut_if_needed(numInputGVCFs, numMergedFields);
    }
    private:
    //why not unordered_map? because I feel the need, the need for speed
    vector<vector<int>> m_inputs_2_merged_lut;
    vector<vector<int>> m_merged_2_inputs_lut;
    void reset_vector(vector<int>& vec, unsigned from=0u);
    void resize_and_reset_vector(vector<int>& vec, unsigned new_size);
    void resize_and_reset_lut(vector<vector<int>>& lut, unsigned new_lut_size, unsigned new_size, unsigned& numRowsVar, unsigned& numColsVar);

    inline int get_lut_value(vector<vector<int>>& lut, int rowIdx, int columnIdx)
    {
      ASSERT(rowIdx >= 0);
      ASSERT(rowIdx < static_cast<int>(lut.size()));
      ASSERT(columnIdx >= 0);
      ASSERT(columnIdx < static_cast<int>(lut[rowIdx].size()));
      return lut[rowIdx][columnIdx];
    }

    inline void set_lut_value(vector<vector<int>>& lut, int rowIdx, int columnIdx, int value)
    {
      ASSERT(rowIdx >= 0);
      ASSERT(rowIdx < static_cast<int>(lut.size()));
      ASSERT(columnIdx >= 0);
      ASSERT(columnIdx < static_cast<int>(lut[rowIdx].size()));
      lut[rowIdx][columnIdx] = value;
    }

    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    inline void set_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    { set_lut_value(m_inputs_2_merged_lut, inputGVCFIdx, inputIdx, mergedIdx); } 

    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    inline void set_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    { set_lut_value(m_inputs_2_merged_lut, inputIdx, inputGVCFIdx, mergedIdx); } 

    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    inline void set_input_idx_for_merged(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    { set_lut_value(m_merged_2_inputs_lut, inputGVCFIdx, mergedIdx, inputIdx); }

    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    inline void set_input_idx_for_merged(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    { set_lut_value(m_merged_2_inputs_lut, mergedIdx, inputGVCFIdx, inputIdx); }

  };

  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  class MergedVCFAllelesIdxLUT
  : public MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>
  {
    public:
      MergedVCFAllelesIdxLUT(unsigned numInputGVCFs);
      void resize_luts_if_needed(unsigned numMergedAlleles);
    private:
      unsigned m_max_num_alleles;
  };

  /*NOTE: Needs explicit instantiation in .cpp file to use this typedef*/
  typedef MergedVCFAllelesIdxLUT<true,true> CombineAllelesLUT;

  //This is the VariantHeaderMerger class (see typedef below)
  //Templated to optimize memory layout
#define DEFAULT_INIT_NUM_INPUT_VCFS 10u
#define DEFAULT_INIT_NUM_FIELDS 30u
#define DEFAULT_INIT_NUM_SAMPLES 10u
  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  class MergedVCFHeaderFieldsIdxLUT
  {
    public:
      MergedVCFHeaderFieldsIdxLUT()
	: m_header_fields_LUT{ DEFAULT_INIT_NUM_INPUT_VCFS,  DEFAULT_INIT_NUM_FIELDS },
	m_samples_LUT { DEFAULT_INIT_NUM_INPUT_VCFS, DEFAULT_INIT_NUM_SAMPLES }
      {
	reset();
	m_num_input_vcfs_allocated = DEFAULT_INIT_NUM_INPUT_VCFS;
	m_num_merged_fields_allocated = DEFAULT_INIT_NUM_FIELDS;
	m_num_merged_samples_allocated = DEFAULT_INIT_NUM_SAMPLES;
      }
      MergedVCFHeaderFieldsIdxLUT(const shared_ptr<bcf_hdr_t>& merged_header,
	  const std::vector<std::shared_ptr<bcf_hdr_t>>& input_vcf_headers)
	: m_header_fields_LUT{ static_cast<unsigned>(input_vcf_headers.size()), static_cast<unsigned>(merged_header->n[BCF_DT_ID]) },
	m_samples_LUT { static_cast<unsigned>(input_vcf_headers.size()), static_cast<unsigned>(bcf_hdr_nsamples(merged_header.get())) }
      {
	reset();
	m_num_input_vcfs_allocated = input_vcf_headers.size();
	m_num_merged_fields_allocated = merged_header->n[BCF_DT_ID];
	m_num_merged_samples_allocated = bcf_hdr_nsamples(merged_header.get());
	m_merged_vcf_header_ptr = merged_header.get();
	auto input_vcf_idx = 0u;
	for(auto input_header : input_vcf_headers)
	{
	  auto input_header_ptr = input_header.get();
	  m_input_vcf_headers.push_back(input_header_ptr);
	  add_header_fields_mapping(input_header_ptr, input_vcf_idx);
	  add_samples_mapping(input_header_ptr, input_vcf_idx);
	  ++input_vcf_idx;
	}
      }
      ~MergedVCFHeaderFieldsIdxLUT() { clear(); }
      void reset()
      {
	m_input_vcf_headers.clear();
	m_sample2idx_merged.clear();
	m_merged_field_idx_enum_vector.clear();
	m_merged_vcf_header_ptr = nullptr;
	m_num_merged_fields_allocated = 0u;
	m_num_merged_samples_allocated = 0u;
	m_num_input_vcfs_allocated = 0u;
      }
      void clear()
      {
	reset();
	m_header_fields_LUT.clear();
	m_samples_LUT.clear();
      }
      inline int get_merged_sample_idx_for_input(unsigned inputGVCFIdx, int inputSampleIdx)
      { return m_samples_LUT.get_merged_idx_for_input(inputGVCFIdx, inputSampleIdx); }
      inline int get_merged_header_idx_for_input(unsigned inputGVCFIdx, int inputIdx)
      { return m_header_fields_LUT.get_merged_idx_for_input(inputGVCFIdx, inputIdx); }

      inline int get_input_sample_idx_for_merged(unsigned inputGVCFIdx, int mergedSampleIdx)
      { return m_samples_LUT.get_input_idx_for_merged(inputGVCFIdx, mergedSampleIdx); }
      inline int get_input_header_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx)
      { return m_header_fields_LUT.get_input_idx_for_merged(inputGVCFIdx, mergedIdx); }
      void store_merged_field_idx_for_enum(const string& field, unsigned field_enum_idx);
      inline int get_merged_field_idx_for_enum(unsigned field_enum_idx)
      {
	ASSERT(field_enum_idx < m_merged_field_idx_enum_vector.size());
	return m_merged_field_idx_enum_vector[field_enum_idx];
      }
      void add_header(shared_ptr<bcf_hdr_t>& merged_hdr, const shared_ptr<bcf_hdr_t>& hdr);
      static inline bool is_missing(int value)
      { return value == MergedVCFLUTBase<true, true>::m_lut_missing_value; }
    private:
      MergedVCFLUTBase<fields_forward_LUT_ordering, fields_reverse_LUT_ordering> m_header_fields_LUT;
      MergedVCFLUTBase<samples_forward_LUT_ordering, samples_reverse_LUT_ordering> m_samples_LUT;
      //Header fields mapping
      void add_header_fields_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx);
      //Samples mapping
      void add_samples_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx);
      //Global sample names to idx mapping
      unordered_map<string,int> m_sample2idx_merged;
      void resize_luts_if_needed();

      vector<bcf_hdr_t*> m_input_vcf_headers;
      vector<int> m_merged_field_idx_enum_vector;
      bcf_hdr_t* m_merged_vcf_header_ptr;
      unsigned m_num_merged_fields_allocated;
      unsigned m_num_merged_samples_allocated;
      unsigned m_num_input_vcfs_allocated;
  };
  /*NOTE: Needs explicit instantiation in .cpp file to use this typedef*/
  typedef MergedVCFHeaderFieldsIdxLUT<true,true,true,true> VariantHeaderMerger;
}

#endif
