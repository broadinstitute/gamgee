set(SOURCE_FILES
    base_quals.cpp
    base_quals.h
    cigar.cpp
    cigar.h
    exceptions.h
    fastq.cpp
    fastq.h
    fastq_iterator.cpp
    fastq_iterator.h
    fastq_reader.cpp
    fastq_reader.h
    gamgee.h
    genotype.cpp
    genotype.h
    indexed_sam_iterator.cpp
    indexed_sam_iterator.h
    indexed_sam_reader.h
    indexed_variant_iterator.cpp
    indexed_variant_iterator.h
    indexed_variant_reader.h
    individual_field.h
    individual_field_iterator.h
    individual_field_value.h
    individual_field_value_iterator.h
    interval.cpp
    interval.h
    missing.h
    multiple_variant_iterator.cpp
    multiple_variant_iterator.h
    multiple_variant_reader.h
    read_bases.cpp
    read_bases.h
    read_group.cpp
    read_group.h
    reference_iterator.cpp
    reference_iterator.h
    reference_map.cpp
    reference_map.h
    sam_builder.cpp
    sam_builder_data_field.cpp
    sam_builder_data_field.h
    sam_builder.h
    sam.cpp
    sam.h
    sam_header.cpp
    sam_header.h
    sam_iterator.cpp
    sam_iterator.h
    sam_pair_iterator.cpp
    sam_pair_iterator.h
    sam_reader.h
    sam_tag.h
    sam_writer.cpp
    sam_writer.h
    shared_field.h
    shared_field_iterator.h
    synced_variant_iterator.cpp
    synced_variant_iterator.h
    synced_variant_reader.h
    utils/file_utils.cpp
    utils/file_utils.h
    utils/genotype_utils.cpp
    utils/genotype_utils.h
    utils/hts_memory.cpp
    utils/hts_memory.h
    utils/short_value_optimized_storage.h
    utils/utils.cpp
    utils/utils.h
    utils/variant_field_type.cpp
    utils/variant_field_type.h
    utils/variant_utils.cpp
    utils/variant_utils.h
    variant_builder.cpp
    variant_builder.h
    variant_builder_individual_field.h
    variant_builder_individual_region.cpp
    variant_builder_individual_region.h
    variant_builder_shared_region.cpp
    variant_builder_shared_region.h
    variant.cpp
    variant_filters.h
    variant_filters_iterator.h
    variant.h
    variant_header_builder.cpp
    variant_header_builder.h
    variant_header.cpp
    variant_header.h
    variant_iterator.cpp
    variant_iterator.h
    variant_reader.h
    variant_writer.cpp
    variant_writer.h
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
