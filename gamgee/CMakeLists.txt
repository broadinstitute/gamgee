set(SOURCE_FILES
    sam/base_quals.cpp
    sam/base_quals.h
    sam/cigar.cpp
    sam/cigar.h
    exceptions.h
    fastq.cpp
    fastq.h
    fastq_iterator.cpp
    fastq_iterator.h
    fastq_reader.cpp
    fastq_reader.h
    gamgee.h
    variant/genotype.cpp
    variant/genotype.h
    sam/indexed_sam_iterator.cpp
    sam/indexed_sam_iterator.h
    sam/indexed_sam_reader.h
    variant/indexed_variant_iterator.cpp
    variant/indexed_variant_iterator.h
    variant/indexed_variant_reader.h
    variant/individual_field.h
    variant/individual_field_iterator.h
    variant/individual_field_value.h
    variant/individual_field_value_iterator.h
    interval.cpp
    interval.h
    missing.h
    variant/multiple_variant_iterator.cpp
    variant/multiple_variant_iterator.h
    variant/multiple_variant_reader.h
    variant/variant_header_merger.h
    variant/variant_header_merger.cpp
    sam/header_comment.cpp
    sam/header_comment.h
    sam/header_line.cpp
    sam/header_line.h
    sam/read_bases.cpp
    sam/read_bases.h
    sam/read_group.cpp
    sam/read_group.h
    sam/program.cpp
    sam/program.h
    variant/reference_block_splitting_variant_iterator.cpp
    variant/reference_block_splitting_variant_iterator.h
    reference_iterator.cpp
    reference_iterator.h
    reference_map.cpp
    reference_map.h
    sam/sam_header_builder.cpp
    sam/sam_header_builder.h
    sam/sam_builder.cpp
    sam/sam_builder_data_field.cpp
    sam/sam_builder_data_field.h
    sam/sam_builder.h
    sam/sam.cpp
    sam/sam.h
    sam/sam_header.cpp
    sam/sam_header.h
    sam/sam_iterator.cpp
    sam/sam_iterator.h
    sam/sam_pair_iterator.cpp
    sam/sam_pair_iterator.h
    sam/sam_reader.h
    sam/sam_tag.cpp
    sam/sam_tag.h
    sam/sam_writer.cpp
    sam/sam_writer.h
    variant/shared_field.h
    variant/shared_field_iterator.h
    variant/synced_variant_iterator.cpp
    variant/synced_variant_iterator.h
    variant/synced_variant_reader.h
    utils/file_utils.cpp
    utils/file_utils.h
    utils/genotype_utils.cpp
    utils/genotype_utils.h
    utils/hts_memory.cpp
    utils/hts_memory.h
    utils/sam_tag_utils.cpp
    utils/sam_tag_utils.h
    utils/short_value_optimized_storage.h
    utils/utils.cpp
    utils/utils.h
    utils/variant_field_type.cpp
    utils/variant_field_type.h
    utils/variant_utils.cpp
    utils/variant_utils.h
    utils/merged_vcf_lut.h
    utils/merged_vcf_lut.cpp
    variant/variant_builder.cpp
    variant/variant_builder.h
    variant/variant_builder_individual_field.h
    variant/variant_builder_individual_region.cpp
    variant/variant_builder_individual_region.h
    variant/variant_builder_multi_sample_vector.h
    variant/variant_builder_shared_region.cpp
    variant/variant_builder_shared_region.h
    variant/variant.cpp
    variant/variant_filters.h
    variant/variant_filters_iterator.h
    variant/variant.h
    variant/variant_header_builder.cpp
    variant/variant_header_builder.h
    variant/variant_header.cpp
    variant/variant_header.h
    variant/variant_iterator.cpp
    variant/variant_iterator.h
    variant/variant_reader.h
    variant/variant_writer.cpp
    variant/variant_writer.h
    zip.h
    )

add_library(gamgee STATIC ${SOURCE_FILES})
add_dependencies(gamgee htslib)

# set install paths for the library binaries
install(TARGETS gamgee
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        COMPONENT library)

# set install paths for headers
# note: this installs all headers; if there are any private headers, we may need to list them explicitly
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
        FILES_MATCHING PATTERN "*.h")
