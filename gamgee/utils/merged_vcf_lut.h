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

  /**
   * LUT = Look Up Table (to avoid confusion with map, unordered_map etc)
   * @brief Base class to store look up information between fields of merged header and input headers
   * Each MergedVCFLUTBase object contains 2 tables (vector of vector): one for mapping input field idx to merged field idx (m_inputs_2_merged_lut)
   * and the second for mapping merged field idx to input field idx (m_merged_2_inputs_lut).
   *
   * Missing field information is stored as -1, but should be checked with is_missing() function
   * 
   * The boolean template parameters specify how the 2 tables are laid out in memory - whether the outer vector corresponds to fields or input vcfs.
   * For example, in object of type MergedVCFLUTBase<true, true>, both LUTs are laid out such that m_inputs_2_merged_lut[0] contains mappings 
   * for all fields for input VCF file 0. This would lead to fast traversal of all fields for a given input VCF (cache locality). 
   * However, traversing over all input VCFs for a given field would be slow (many cache misses).
   * The object MergedVCFLUTBase<false,false> would have the exact opposite behavior
   * 
   * The 'best' value of the template parameters depends on the application using the LUT.
   * Almost all the 'complexity' of the code comes from being able to handle the different layouts in a transparent manner
   *
   * Although the class provides functions to resize the tables, for obtaining good performance, reallocations should be extremely
   * infrequent. Making the resize_luts_if_needed() a protected member forces developers to think twice instead of blindly calling this function.
   *
   * Uses the enable_if trick from http://en.cppreference.com/w/cpp/types/enable_if  (foo3) to handle different memory layouts
   * 
   **/
  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  class MergedVCFLUTBase
  {
    template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
    friend class MergedVCFHeaderFieldsIdxLUT;
    public:
    /**
     * @brief: clear all mappings
     */
    inline void reset_luts()
    {
      for(auto& vec : m_inputs_2_merged_lut)
	reset_vector(vec);
      for(auto& vec : m_merged_2_inputs_lut)
	reset_vector(vec);
    }

    /*
     * @brief Add a valid mapping between input VCF and merged VCF
     * @note all parameters should be valid parameters, no -1 allowed, use reset_() functions to invalidate existing mapping
     * @param inputGVCFIdx index of the input VCF file
     * @param inputIdx index of the field in the input VCF file - field could be anything header field,sample,allele etc
     * @param mergedIdx index of the field in the merged VCF file
     */
    inline void add_input_merged_idx_pair(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    {
      set_merged_idx_for_input(inputGVCFIdx, inputIdx, mergedIdx);
      set_input_idx_for_merged(inputGVCFIdx, inputIdx, mergedIdx);
    }

    /**
     * @brief Get field idx for input VCF inputGVCFIdx corresponding to field idx mergedIdx in the mergedVCF file
     * @note Uses the enable_if trick from http://en.cppreference.com/w/cpp/types/enable_if  (foo3) to handle different memory layouts
     * The enable_if<M> corresponds to the case where merged_2_inputs_LUT_is_input_ordered = true, hence, the rows correspond to input VCFs
     * The enable_if<!M> corresponds to the case where merged_2_inputs_LUT_is_input_ordered = false, hence, the rows correspond to fields
     * @param inputGVCFIdx index of the input VCF file
     * @param mergedIdx index of the field in the merged VCF file
     * @return index of the field in the input VCF file
     */
    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    inline int get_input_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx)
    { return get_lut_value(m_merged_2_inputs_lut, inputGVCFIdx, mergedIdx); }
    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    inline int get_input_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx)
    { return get_lut_value(m_merged_2_inputs_lut, mergedIdx, inputGVCFIdx); }

    /**
     * @brief Get field idx for the merged VCF corresponding to field idx inputIdx in the input VCF of index inputGVCFIdx
     * @note Uses the enable_if trick from http://en.cppreference.com/w/cpp/types/enable_if  (foo3) to handle different memory layouts
     * The enable_if<M> corresponds to the case where inputs_2_merged_LUT_is_input_ordered = true, hence, the rows correspond to input VCFs
     * The enable_if<!M> corresponds to the case where inputs_2_merged_LUT_is_input_ordered = false, hence, the rows correspond to fields
     * @param inputGVCFIdx index of the input VCF file
     * @param inputIdx index of the field in the input VCF file
     * @return index of the field in the merged VCF file
     */
    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    inline int get_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx)
    { return get_lut_value(m_inputs_2_merged_lut, inputGVCFIdx, inputIdx); }
    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    inline int get_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx)
    { return get_lut_value(m_inputs_2_merged_lut, inputIdx, inputGVCFIdx); }

    /**
     * @brief reset/invalidate merged field index for field inputIdx of input VCF inputGVCFIdx
     * @param inputGVCFIdx index of the input VCF file
     * @param inputIdx index of the field in the input VCF file
     */
    inline void reset_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx)
    {	set_merged_idx_for_input(inputGVCFIdx, inputIdx, MergedVCFLUTBase::m_lut_missing_value); }
    /**
     * @brief reset/invalidate the input field index for input VCF inputGVCFIdx for merged field mergedIdx
     * @param inputGVCFIdx index of the input VCF file
     * @param mergedIdx index of the field in the merged VCF file
     */
    inline void reset_input_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx)
    {	set_input_idx_for_merged(inputGVCFIdx, MergedVCFLUTBase::m_lut_missing_value, mergedIdx); }

    /**
     * @brief check if the value returned by get_*() functions for the LUT is invalid
     * @param lut_val value obtained by get_*() functions of the LUT
     * @return true if lut_val is invalid, else false
     */
    static inline bool is_missing(const int lut_val) {  return (lut_val == MergedVCFLUTBase::m_lut_missing_value); }
    protected:
    static const int m_lut_missing_value = -1;
    //Only inherited classes should call constructor,destructor etc
    MergedVCFLUTBase(); 
    MergedVCFLUTBase(unsigned numInputGVCFs, unsigned numMergedFields);
    ~MergedVCFLUTBase() { clear(); }
    /**
     * @brief deallocates memory
     */
    void clear();

    unsigned m_num_input_vcfs;
    unsigned m_num_merged_fields;

    /**
     *  @brief resize LUT functions 
     *  @note should be called relatively infrequently (more precisely, the reallocation code inside these resize functions should be called
     *  infrequently
     *  @param numInputGVCFs number of input VCFs
     *  @param numMergedFields number of fields combined across all input VCFs
     */
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

    /*
     * @brief wrapper around single LUT resize functions
     */
    void resize_luts_if_needed(unsigned numInputGVCFs, unsigned numMergedFields)
    {
      resize_merged_2_inputs_lut_if_needed(numInputGVCFs, numMergedFields);
      resize_inputs_2_merged_lut_if_needed(numInputGVCFs, numMergedFields);
    }
    private:
    //why not unordered_map? because I feel the need, the need for speed
    vector<vector<int>> m_inputs_2_merged_lut;
    vector<vector<int>> m_merged_2_inputs_lut;
    /**
     * @brief invalidate/reset all mappings in a vector
     * @note sets all elements to -1
     * @param vec the vector to reset
     * @param from offset in the vector from which to start reset, 0 by default
     */
    void reset_vector(vector<int>& vec, unsigned from=0u);
    /**
     * @brief resize and reset a vector
     * @note resize and reset is done only if new_size > vec.size()
     */
    void resize_and_reset_vector(vector<int>& vec, unsigned new_size);
    /**
     * @brief resize and reset a LUT
     * @note resize and reset is done only if new_size > old_size
     */
    void resize_and_reset_lut(vector<vector<int>>& lut, unsigned new_lut_size, unsigned new_size, unsigned& numRowsVar, unsigned& numColsVar);

    /**
     * @brief get LUT value at a particular row,column
     * @note should be called only from the public wrapper functions get_*() as the wrappers take care of memory layout
     * @param lut LUT to access
     * @param rowIdx row
     * @param columnIdx column
     * @return value at lut[row][column], could be invalid, check with is_missing()
     */
    inline int get_lut_value(vector<vector<int>>& lut, int rowIdx, int columnIdx)
    {
      ASSERT(rowIdx >= 0);
      ASSERT(rowIdx < static_cast<int>(lut.size()));
      ASSERT(columnIdx >= 0);
      ASSERT(columnIdx < static_cast<int>(lut[rowIdx].size()));
      return lut[rowIdx][columnIdx];
    }
    
    /**
     * @brief set LUT value at a particular row,column
     * @note should be called only from the public wrapper functions add_input_merged_idx_pair() or reset_*() as the wrappers take care of memory layout
     * @param lut LUT to access
     * @param rowIdx row
     * @param columnIdx column
     * @param value value to write at lut[row][column] 
     */
    inline void set_lut_value(vector<vector<int>>& lut, int rowIdx, int columnIdx, int value)
    {
      ASSERT(rowIdx >= 0);
      ASSERT(rowIdx < static_cast<int>(lut.size()));
      ASSERT(columnIdx >= 0);
      ASSERT(columnIdx < static_cast<int>(lut[rowIdx].size()));
      lut[rowIdx][columnIdx] = value;
    }

    /**
     * @brief set merged field idx value (mergedIdx) corresponding to field idx inputIdx for input VCF inputGVCFIdx
     * @note should be called only from the public wrapper functions add_input_merged_idx_pair() or reset_*() as the wrappers take care of memory layout
     * @param inputGVCFIdx index of the input VCF file
     * @param inputIdx index of the field in the input VCF file
     * @param mergedIdx index of the field in the merged VCF file
     */
    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    inline void set_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    { set_lut_value(m_inputs_2_merged_lut, inputGVCFIdx, inputIdx, mergedIdx); } 

    template <bool M = inputs_2_merged_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    inline void set_merged_idx_for_input(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    { set_lut_value(m_inputs_2_merged_lut, inputIdx, inputGVCFIdx, mergedIdx); } 

    /**
     * @brief set input field idx value (inputIdx) for input VCF inputGVCFIdx corresponding to field idx mergedIdx in the merged VCF
     * @note should be called only from the public wrapper functions add_input_merged_idx_pair() or reset_*() as the wrappers take care of memory layout
     * @param inputGVCFIdx index of the input VCF file
     * @param inputIdx index of the field in the input VCF file
     * @param mergedIdx index of the field in the merged VCF file
     */
    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<M>::type* = nullptr>
    inline void set_input_idx_for_merged(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    { set_lut_value(m_merged_2_inputs_lut, inputGVCFIdx, mergedIdx, inputIdx); }

    template <bool M = merged_2_inputs_LUT_is_input_ordered, typename std::enable_if<!M>::type* = nullptr>
    inline void set_input_idx_for_merged(unsigned inputGVCFIdx, int inputIdx, int mergedIdx)
    { set_lut_value(m_merged_2_inputs_lut, mergedIdx, inputGVCFIdx, inputIdx); }

  };

  /**
   * @brief LUT class for storing mappings between allele vectors in the merged file and input VCF files
   * Since the #alleles per site is expected to be small, this class sets the number of fields to 10. This makes any subsequent re-allocations
   * unlikely. The function resize_luts_if_needed() will almost always return immediately after failing the if condition
   */
  template<bool inputs_2_merged_LUT_is_input_ordered, bool merged_2_inputs_LUT_is_input_ordered>
  class MergedVCFAllelesIdxLUT
  : public MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>
  {
#define DEFAULT_INIT_NUM_ALLELES 10u
    public:
      MergedVCFAllelesIdxLUT(unsigned numInputGVCFs)
	: MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>(numInputGVCFs, DEFAULT_INIT_NUM_ALLELES)
      { m_max_num_alleles = DEFAULT_INIT_NUM_ALLELES; }
      inline void resize_luts_if_needed(unsigned numMergedAlleles)
      {
	if(numMergedAlleles > m_max_num_alleles)
	{
	  MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::resize_luts_if_needed(
	      MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>::m_num_input_vcfs, numMergedAlleles); 
	  m_max_num_alleles = numMergedAlleles;
	}
      }
    private:
      unsigned m_max_num_alleles;
  };

  /*NOTE: Needs explicit instantiation in .cpp file to use this typedef*/
  typedef MergedVCFAllelesIdxLUT<true,true> CombineAllelesLUT;

  /*
   * VariantHeaderMerger class (see typedef below)
   * @brief Class which maintains mappings between input VCFs and merged VCF header fields and samples. It is templated to select
   * the 'best' memory layout (see notes on MergedVCFLUTBase)
   * The class contains two LUTs - m_header_fields_LUT and m_samples_LUT of type MergedVCFLUTBase<> for storing mapping for
   * header fields (FMT, FLT, INFO) and samples. Additionally, this class provides wrapper functions to extract information
   * from the LUTs
   * This class provides functions to add new fields as well as add new VCF files. Internally, it calls resize functions as and when needed
   */
#define DEFAULT_INIT_NUM_INPUT_VCFS 10u
#define DEFAULT_INIT_NUM_FIELDS 30u
#define DEFAULT_INIT_NUM_SAMPLES 10u
  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  class MergedVCFHeaderFieldsIdxLUT
  {
    public:
      /**
       * @brief empty constructor, initialize 'large' LUTs
       */
      MergedVCFHeaderFieldsIdxLUT()
	: m_header_fields_LUT{ DEFAULT_INIT_NUM_INPUT_VCFS,  DEFAULT_INIT_NUM_FIELDS },
	m_samples_LUT { DEFAULT_INIT_NUM_INPUT_VCFS, DEFAULT_INIT_NUM_SAMPLES }
      {
	reset();
	m_num_input_vcfs_allocated = DEFAULT_INIT_NUM_INPUT_VCFS;
	m_num_merged_fields_allocated = DEFAULT_INIT_NUM_FIELDS;
	m_num_merged_samples_allocated = DEFAULT_INIT_NUM_SAMPLES;
      }
      /**
       * @brief Constructor where the merged VCF's header is already built
       * @param merged_header fully built header of the merged VCF
       * @param input_vcf_headers vector of headers of input VCFs which are being merged
       */
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
      /**
       * @brief: resets all mappings, but does not de-allocate LUT memory
       */
      void reset()
      {
	m_input_vcf_headers.clear();
	m_sample2idx_merged.clear();
	m_merged_field_idx_enum_vector.clear();
	m_merged_vcf_header_ptr = nullptr;
	m_num_merged_fields_allocated = 0u;
	m_num_merged_samples_allocated = 0u;
	m_num_input_vcfs_allocated = 0u;
	m_header_fields_LUT.reset_luts();
	m_samples_LUT.reset_luts();
      }
      /**
       * @brief: resets all mappings, de-allocates LUT memory
       */
      void clear()
      {
	reset();
	m_header_fields_LUT.clear();
	m_samples_LUT.clear();
      }
      /**
       * @brief Get sample idx for the merged VCF corresponding to sample idx inputSampleIdx in the input VCF of index inputGVCFIdx
       * @param inputGVCFIdx index of the input VCF file
       * @param inputSampleIdx index of the sample in the input VCF file
       * @return index of the sample in the merged VCF file
       */
      inline int get_merged_sample_idx_for_input(unsigned inputGVCFIdx, int inputSampleIdx)
      { return m_samples_LUT.get_merged_idx_for_input(inputGVCFIdx, inputSampleIdx); }
      /**
       * @brief Get header field (FLT/FMT/INFO) idx for the merged VCF corresponding to field idx inputIdx in the input VCF of index inputGVCFIdx
       * @param inputGVCFIdx index of the input VCF file
       * @param inputIdx index of the field in the input VCF file
       * @return index of the field in the merged VCF file
       */
      inline int get_merged_header_idx_for_input(unsigned inputGVCFIdx, int inputIdx)
      { return m_header_fields_LUT.get_merged_idx_for_input(inputGVCFIdx, inputIdx); }
      /**
       * @brief Get sample idx for the input VCF inputGVCFIdx corresponding to sample mergedSampleIdx in the merged VCF 
       * @param inputGVCFIdx index of the input VCF file
       * @param mergedSampleIdx index of the sample in the merged VCF file
       * @return index of the sample in the input VCF file inputGVCFIdx
       */
      inline int get_input_sample_idx_for_merged(unsigned inputGVCFIdx, int mergedSampleIdx)
      { return m_samples_LUT.get_input_idx_for_merged(inputGVCFIdx, mergedSampleIdx); }
      /**
       * @brief Get header field (FLT/FMT/INFO) idx for the input VCF inputGVCFIdx corresponding to field mergedIdx in the merged VCF 
       * @param inputGVCFIdx index of the input VCF file
       * @param mergedIdx index of the field in the merged VCF file
       * @return index of the field in the input VCF file inputGVCFIdx
       */
      inline int get_input_header_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx)
      { return m_header_fields_LUT.get_input_idx_for_merged(inputGVCFIdx, mergedIdx); }
      /**
       * @brief add a new header into the merged header and update LUTs
       * @param merged_hdr  current merged header
       * @param hdr new input header to add
       */
      void add_header(shared_ptr<bcf_hdr_t>& merged_hdr, const shared_ptr<bcf_hdr_t>& hdr);
      /**
       * @brief check if the value returned by get_*() functions for the LUT is invalid
       * @param lut_val value obtained by get_*() functions of the LUT
       * @return true if lut_val is invalid, else false
       */
      static inline bool is_missing(int value)
      { return value == MergedVCFLUTBase<true, true>::m_lut_missing_value; }

      /**
       * @brief utility function for storing index of frequently used fields in the merged VCF
       * Sometimes the user/developer may know beforehand that certain fields are needed/accesssed (for example "PL")
       * for every variant in a merged VCF
       * Instead of doing a string search everytime, the user may wish to define an enum corresponding to the relevant fields
       * and store the index of the relevant fields in the merged VCF after the merged header is built. Subsequent accesses to
       * the fields need not use any string searches, instead directly using the idx in the LUT
       * Example, the user could define an enum { PL_FIELD=0, AD_FIELD, AF_FIELD } and call the function
       * store_merged_field_idx_for_enum("PL", PL_FIELD);
       * store_merged_field_idx_for_enum("AD", AD_FIELD); ...
       *
       * While processing variants, the user could directly use the merged field idx using the function
       * auto merged_PL_idx = get_merged_field_idx_for_enum(PL_FIELD);
       * This completely avoids the string searches during variant processing
       */
      void store_merged_field_idx_for_enum(const string& field, unsigned field_enum_idx);
      inline int get_merged_field_idx_for_enum(unsigned field_enum_idx)
      {
	ASSERT(field_enum_idx < m_merged_field_idx_enum_vector.size());
	return m_merged_field_idx_enum_vector[field_enum_idx];
      }
    private:
      /**
       * @brief function to resize LUTs if needed
       */
      void resize_luts_if_needed();
      //LUT for VCF header fields (FMT/FLT/INFO)
      MergedVCFLUTBase<fields_forward_LUT_ordering, fields_reverse_LUT_ordering> m_header_fields_LUT;
      //LUT for samples
      MergedVCFLUTBase<samples_forward_LUT_ordering, samples_reverse_LUT_ordering> m_samples_LUT;
      //Header fields mapping
      void add_header_fields_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx);
      //Samples mapping
      void add_samples_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx);
      //Global sample names to idx mapping
      unordered_map<string,int> m_sample2idx_merged;
      //Input VCF headers
      vector<bcf_hdr_t*> m_input_vcf_headers;
      //Merged header
      bcf_hdr_t* m_merged_vcf_header_ptr;
      //sizes of the LUTs - to determine when to reallocate
      unsigned m_num_merged_fields_allocated;
      unsigned m_num_merged_samples_allocated;
      unsigned m_num_input_vcfs_allocated;
      //vector to store merged field idxs corresponding to user-defined enum
      vector<int> m_merged_field_idx_enum_vector;
  };
  /*NOTE: Needs explicit instantiation in .cpp file to use this typedef*/
  typedef MergedVCFHeaderFieldsIdxLUT<true,true,true,true> VariantHeaderMerger;
}

#endif
