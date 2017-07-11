set(SOURCE_FILES
    cigar_test.cpp
    fastq_reader_test.cpp
    fastq_test.cpp
    genotypes_test.cpp
    indexed_sam_reader_test.cpp
    indexed_variant_reader_test.cpp
    interval_test.cpp
    main.cpp
    missing_test.cpp
    multiple_variant_reader_test.cpp
    header_comment_test.cpp
    header_line_test.cpp
    read_group_test.cpp
    program_test.cpp
    reference_block_splitting_variant_reader_test.cpp
    reference_test.cpp
    sam_header_builder_test.cpp
    sam_builder_test.cpp
    sam_header_test.cpp
    sam_reader_test.cpp
    sam_test.cpp
    select_if_test.cpp
    short_value_optimized_storage_test.cpp
    synced_variant_reader_test.cpp
    test_utils.h
    utils_test.cpp
    variant_builder_multi_sample_vector_test.cpp
    variant_builder_test.cpp
    variant_header_test.cpp
    variant_reader_test.cpp
    variant_test.cpp)

add_executable(gamgee_test EXCLUDE_FROM_ALL ${SOURCE_FILES})

target_compile_definitions(gamgee_test PUBLIC -DBOOST_TEST_DYN_LINK)
target_link_libraries(gamgee_test gamgee ${htslib_LIB} ${Boost_LIBRARIES} pthread z bz2 lzma)
add_dependencies(gamgee_test htslib)

add_custom_target(run_test COMMAND ${CMAKE_BINARY_DIR}/test/gamgee_test DEPENDS gamgee_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
