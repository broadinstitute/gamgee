#ifndef __gamgee_variant_header_merger__
#define __gamgee_variant_header_merger__

#include <boost/shared_ptr.hpp>
#include <string>
#include <unordered_map>
#include "utils/merged_vcf_lut.h"

namespace gamgee
{
  /*
   * VariantHeaderMerger class
   * @brief Class for merging VCF headers and maintain field and sample mappings
   * Class which:
   * (a) Merges multiple input VCF headers to obtain a merged VCF header
   * (b) Creates and updates mappings between input VCFs and merged VCF header fields and samples.
   * 
   * Provides functions to add new VCF headers for merging
   * Example usage:
   * shared_ptr<bcf_hdr_t> vcf_hdr;
   * vector<shared_ptr<bcf_hdr_t>> hdr_vec1;
   * vector<shared_ptr<bcf_hdr_t>> hdr_vec2;
   *
   * Start with empty object, add headers later
   * VariantHeaderMerger X;
   * X.add_header(vcf_hdr);
   * X.add_headers(hdr_vec1);
   * const shared_ptr<bcf_hdr_t>& merged_hdr = X.get_merged_header();
   *
   * Start with 1, add more later
   * VariantHeaderMerger Y { vcf_hdr };
   * Y.add_headers(hdr_vec1);
   *
   * Start with many, add more later
   * VariantHeaderMerger Z { hdr_vec1 };
   * Z.add_headers(hdr_vec2);
   *
   * The class contains two LUTs - m_header_fields_LUT and m_samples_LUT of type MergedVCFLUTBase<> for storing mapping for
   * header fields (FMT, FLT, INFO) and samples respectively. The class is templated to select the 'best' memory layout.
   * 
   * Each of the two MergedVCFLUTBase<> objects (m_header_fields_LUT, m_samples_LUT) contains two matrices (vector<vector<int>>):
   * m_inputs_2_merged_lut and m_merged_2_inputs_lut. The first stores the mapping from input VCF fields to the merged VCF fields while
   * the second stores the mapping in the opposite direction.
   * You can layout each matrix in one of the 2 following ways:
   * (a) matrix[i][j]  corresponds to input VCF i and field j 
   * (b) matrix[i][j]  corresponds to field i and input VCF j
   * Option (a) is optimal where you are looking at all the fields of a VCF in quick succession,
   * while (b) is optimal when you are looking at all VCFs for a particular field.
   * The 2 boolean template parameters for MergedVCFLUTBase<inputs_2_merged_LUT_is_input_ordered, merged_2_inputs_LUT_is_input_ordered>
   * control the layout of the two matrices. If the parameter inputs_2_merged_LUT_is_input_ordered is true, then layout (a) is selected for
   * m_inputs_2_merged_lut, else layout (b)
   *
   * VariantHeaderMerger has 4 boolean parameters - two for the m_header_fields_LUT (fields_*_LUT_ordering) and two for the
   * m_samples_LUT (samples_*_LUT_ordering)
   * fields_forward_LUT_ordering: controls the layout of the input VCF fields to merged VCF fields matrix in m_header_fields_LUT
   * fields_reverse_LUT_ordering: controls the layout of the merged VCF fields to input VCF fields matrix in m_header_fields_LUT
   * samples_forward_LUT_ordering: controls the layout of the input VCF fields to merged VCF samples matrix in m_samples_LUT
   * samples_reverse_LUT_ordering: controls the layout of the merged VCF samples to input VCF samples matrix in m_samples_LUT
   *
   * We create type aliases InputOrderedVariantHeaderMerger and FieldOrderedVariantHeaderMerger for specialized template
   * instantiations of VariantHeaderMerger.
   * (a) InputOrderedVariantHeaderMerger - fast traversal of all field lookups for a single input VCF
   * (b) FieldOrderedVariantHeaderMerger - fast traversal of all input VCFs for a single field lookup
   *
   * Example usage:
   * vector<VariantHeader> hdr_vec;
   * VariantHeaderMerger Z { hdr_vec };
   * auto input_vcf_idx = 0u;
   * auto input_PL_idx = hdr_vec[input_vcf_idx].field_index("PL", BCF_HL_FMT);
   * auto merged_PL_idx = Z.get_merged_header_idx_for_input(input_vcf_idx, input_PL_idx);
   * auto refind_input_PL_idx = Z.get_input_header_idx_for_merged(input_vcf_idx, merged_PL_idx);
   * assert(refind_input_PL_idx == input_PL_idx);
   */

  template<bool fields_forward_LUT_ordering, bool fields_reverse_LUT_ordering, bool samples_forward_LUT_ordering, bool samples_reverse_LUT_ordering>
  class VariantHeaderMerger
  {
    private:
      static const auto m_DEFAULT_INIT_NUM_INPUT_VCFS  = 10u;
      static const auto m_DEFAULT_INIT_NUM_FIELDS = 30u;
      static const auto m_DEFAULT_INIT_NUM_SAMPLES = 10u;
    public:
      /**
       * @brief empty constructor, initialize 'large' LUTs
       */
      VariantHeaderMerger()
        : m_header_fields_LUT{ m_DEFAULT_INIT_NUM_INPUT_VCFS,  m_DEFAULT_INIT_NUM_FIELDS },
        m_samples_LUT { m_DEFAULT_INIT_NUM_INPUT_VCFS, m_DEFAULT_INIT_NUM_SAMPLES },
        m_merged_field_idx_enum_lut { 1u, m_DEFAULT_INIT_NUM_FIELDS }
      {
        reset();
        m_num_input_vcfs_allocated = m_DEFAULT_INIT_NUM_INPUT_VCFS;
        m_num_merged_fields_allocated = m_DEFAULT_INIT_NUM_FIELDS;
        m_num_merged_samples_allocated = m_DEFAULT_INIT_NUM_SAMPLES;
        m_num_enums_allocated = m_DEFAULT_INIT_NUM_FIELDS;
      }
      /**
       * @brief Constructor with a single input VCF header as input 
       * @param input_vcf_header header of input VCF
       */
      VariantHeaderMerger(const std::shared_ptr<bcf_hdr_t>& input_vcf_header)
        : VariantHeaderMerger() 
      {
        add_header(input_vcf_header);
      }
      /**
       * @brief Constructor with a single input VCF header as input 
       * @param input_vcf_header header of input VCF
       */
      VariantHeaderMerger(const VariantHeader& input_vcf_header)
        : VariantHeaderMerger(input_vcf_header.m_header) 
      { }
      /**
       * @brief Constructor with a vector of input VCF headers to be merged
       * @param input_vcf_headers vector of headers of input VCFs which are being merged
       */
      VariantHeaderMerger(const std::vector<std::shared_ptr<bcf_hdr_t>>& input_vcf_headers)
        : m_header_fields_LUT{ static_cast<unsigned>(input_vcf_headers.size()),  m_DEFAULT_INIT_NUM_FIELDS },
        m_samples_LUT { static_cast<unsigned>(input_vcf_headers.size()), m_DEFAULT_INIT_NUM_SAMPLES },
        m_merged_field_idx_enum_lut { 1u, m_DEFAULT_INIT_NUM_FIELDS }
      {
        reset();
        m_num_input_vcfs_allocated = input_vcf_headers.size();
        m_num_merged_fields_allocated = m_DEFAULT_INIT_NUM_FIELDS;
        m_num_merged_samples_allocated = m_DEFAULT_INIT_NUM_SAMPLES;
        m_num_enums_allocated = m_DEFAULT_INIT_NUM_FIELDS;
        add_headers(input_vcf_headers);
      }
      /**
       * @brief Constructor with a vector of input VCF headers to be merged
       * @param input_vcf_headers vector of headers of input VCFs which are being merged
       */
      VariantHeaderMerger(const std::vector<VariantHeader>& input_vcf_headers)
        : m_header_fields_LUT{ static_cast<unsigned>(input_vcf_headers.size()),  m_DEFAULT_INIT_NUM_FIELDS },
        m_samples_LUT { static_cast<unsigned>(input_vcf_headers.size()), m_DEFAULT_INIT_NUM_SAMPLES },
        m_merged_field_idx_enum_lut { 1u, m_DEFAULT_INIT_NUM_FIELDS }
      {
        reset();
        m_num_input_vcfs_allocated = input_vcf_headers.size();
        m_num_merged_fields_allocated = m_DEFAULT_INIT_NUM_FIELDS;
        m_num_merged_samples_allocated = m_DEFAULT_INIT_NUM_SAMPLES;
        m_num_enums_allocated = m_DEFAULT_INIT_NUM_FIELDS;
        add_headers(input_vcf_headers);
      }
      /*
       * @brief No anticipated use of a deep copy for VariantHeaderMerger
       */
      VariantHeaderMerger(const VariantHeaderMerger&) = delete;
      VariantHeaderMerger& operator=(const VariantHeaderMerger&) = delete;
      /*
       * @brief default move constructor for VariantHeaderMerger
       */
      VariantHeaderMerger(VariantHeaderMerger&&) = default;
      VariantHeaderMerger& operator=(VariantHeaderMerger&&) = default;

      ~VariantHeaderMerger() = default;
      /**
       * @brief: resets all mappings, but does not de-allocate LUT memory
       */
      void reset()
      {
        m_input_vcf_headers.clear();
        m_sample2idx_merged.clear();
        m_merged_vcf_header_ptr = nullptr;
        m_num_merged_fields_allocated = 0u;
        m_num_merged_samples_allocated = 0u;
        m_num_input_vcfs_allocated = 0u;
        m_num_enums_allocated = 0u;
        m_header_fields_LUT.reset_luts();
        m_samples_LUT.reset_luts();
        m_merged_field_idx_enum_lut.reset_luts();
      }
      /**
       * @brief: resets all mappings, de-allocates LUT memory
       */
      void clear()
      {
        reset();
        m_header_fields_LUT.clear();
        m_samples_LUT.clear();
        m_merged_field_idx_enum_lut.clear();
      }
      /**
       * @brief add a new header into the merged header and update LUTs
       * @param hdr new input header to add
       */
      void add_header(const std::shared_ptr<bcf_hdr_t>& hdr);
      /**
       * @brief add a new header into the merged header and update LUTs
       * @param hdr new input header to add
       */
      void add_header(const VariantHeader& hdr);
      /**
       * @brief add a vector of new VCF headers into the merged header and update LUTs
       * @param headers vector of new input VCF headers to add
       */
      void add_headers(const std::vector<std::shared_ptr<bcf_hdr_t>>& headers);
      /**
       * @brief add a vector of new VCF headers into the merged header and update LUTs
       * @param headers vector of new input VCF headers to add
       */
      void add_headers(const std::vector<VariantHeader>& headers); 
      /**
       * @brief Get merged VCF header shared_ptr
       * @return return the merged VCF header shared_ptr
       */
      const std::shared_ptr<bcf_hdr_t>& get_raw_merged_header() const { return m_merged_vcf_header_ptr; }
      /**
       * @brief Get merged VCF header
       * @return return the merged VCF header
       */
      const VariantHeader get_merged_header() const { return VariantHeader{ m_merged_vcf_header_ptr }; }
      /*LUT functions*/
      /**
       * @brief Get sample idx for the merged VCF corresponding to sample idx inputSampleIdx in the input VCF of index inputGVCFIdx
       * @param inputGVCFIdx index of the input VCF file
       * @param inputSampleIdx index of the sample in the input VCF file
       * @return index of the sample in the merged VCF file
       */
      inline int get_merged_sample_idx_for_input(unsigned inputGVCFIdx, int inputSampleIdx) const
      { return m_samples_LUT.get_merged_idx_for_input(inputGVCFIdx, inputSampleIdx); }
      /**
       * @brief Get header field (FLT/FMT/INFO) idx for the merged VCF corresponding to field idx inputIdx in the input VCF of index inputGVCFIdx
       * @param inputGVCFIdx index of the input VCF file
       * @param inputIdx index of the field in the input VCF file
       * @return index of the field in the merged VCF file
       */
      inline int get_merged_header_idx_for_input(unsigned inputGVCFIdx, int inputIdx) const
      { return m_header_fields_LUT.get_merged_idx_for_input(inputGVCFIdx, inputIdx); }
      /**
       * @brief Get sample idx for the input VCF inputGVCFIdx corresponding to sample mergedSampleIdx in the merged VCF 
       * @param inputGVCFIdx index of the input VCF file
       * @param mergedSampleIdx index of the sample in the merged VCF file
       * @return index of the sample in the input VCF file inputGVCFIdx
       */
      inline int get_input_sample_idx_for_merged(unsigned inputGVCFIdx, int mergedSampleIdx) const
      { return m_samples_LUT.get_input_idx_for_merged(inputGVCFIdx, mergedSampleIdx); }
      /**
       * @brief Get header field (FLT/FMT/INFO) idx for the input VCF inputGVCFIdx corresponding to field mergedIdx in the merged VCF 
       * @param inputGVCFIdx index of the input VCF file
       * @param mergedIdx index of the field in the merged VCF file
       * @return index of the field in the input VCF file inputGVCFIdx
       */
      inline int get_input_header_idx_for_merged(unsigned inputGVCFIdx, int mergedIdx) const
      { return m_header_fields_LUT.get_input_idx_for_merged(inputGVCFIdx, mergedIdx); }

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
       * This completely avoids string searches during variant processing
       *
       * The reverse function (get_enum_for_merged_field_idx) is also useful:
       * for(i=0;i<bcf1_t->n_info;++i)
       * {
       *   switch(VariantHeaderMerger.get_enum_for_merged_field_idx(bcf1_t->d.info[i].key))
       *   {
       *      case PL_FIELD:
       *        do_something_with_PL();
       *        break;
       *      ....
       *   }
       * }
       */
      void store_merged_field_idx_for_enum(const std::string& field, unsigned field_enum_idx);
      inline int get_merged_field_idx_for_enum(unsigned field_enum_idx) const
      {
        return m_merged_field_idx_enum_lut.get_merged_idx_for_input(0u, field_enum_idx);
      }
      inline int get_enum_for_merged_field_idx(int merged_field_idx) const
      {
        return m_merged_field_idx_enum_lut.get_input_idx_for_merged(0u, merged_field_idx);
      }
    private:
      /**
       * @brief function to resize LUTs if needed
       */
      void resize_luts_if_needed();
      //LUT for VCF header fields (FMT/FLT/INFO)
      utils::MergedVCFLUTBase<fields_forward_LUT_ordering, fields_reverse_LUT_ordering> m_header_fields_LUT;
      //LUT for samples
      utils::MergedVCFLUTBase<samples_forward_LUT_ordering, samples_reverse_LUT_ordering> m_samples_LUT;
      //Header fields mapping
      void add_header_fields_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx);
      //Samples mapping
      void add_samples_mapping(bcf_hdr_t* curr_header, unsigned input_vcf_idx);
      //Global sample names to idx mapping
      std::unordered_map<std::string,int> m_sample2idx_merged;
      //Input VCF headers
      std::vector<std::shared_ptr<bcf_hdr_t>> m_input_vcf_headers;
      //Merged header
      std::shared_ptr<bcf_hdr_t> m_merged_vcf_header_ptr;
      //sizes of the LUTs - to determine when to reallocate
      unsigned m_num_merged_fields_allocated;
      unsigned m_num_merged_samples_allocated;
      unsigned m_num_input_vcfs_allocated;
      unsigned m_num_enums_allocated;
      //LUT to store merged field idxs to user-defined enum mappings
      //The two matrices in this LUT are single row matrices - the matrix m_inputs_2_merged_lut will
      //store the mapping from the user defined enum to the merged VCF header field idx, while the
      //matrix m_merged_2_inputs will store the enum corresponding to the merged field
      utils::MergedVCFLUTBase<true, true> m_merged_field_idx_enum_lut;
  };
  /*NOTE: Needs explicit instantiation in .cpp file to use this type alias*/
  using InputOrderedVariantHeaderMerger = VariantHeaderMerger<true,true,true,true>;
  using FieldOrderedVariantHeaderMerger = VariantHeaderMerger<false,false,false,false>;
}

#endif
