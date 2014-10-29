set(SOURCE_FILES
    cigar_test.cpp
    fastq_reader_test.cpp
    fastq_test.cpp
    genotypes_test.cpp
    indexed_sam_reader_test.cpp
    interval_test.cpp
    main.cpp
    missing_test.cpp
    read_group_test.cpp
    reference_test.cpp
    sam_builder_test.cpp
    sam_header_test.cpp
    sam_reader_test.cpp
    sam_test.cpp
    select_if_test.cpp
    short_value_optimized_storage_test.cpp
    utils_test.cpp
    variant_builder_test.cpp
    variant_header_test.cpp
    variant_reader_test.cpp
    variant_test.cpp)

include("../contrib/htslib.cmake")

include_directories(${Boost_INCLUDE_DIRS})

add_executable(gamgee_test ${SOURCE_FILES})
target_compile_definitions(gamgee_test PUBLIC -DBOOST_TEST_DYN_LINK)
target_link_libraries(gamgee_test gamgee ${htslib_LIB} ${Boost_LIBRARIES} pthread z)
add_dependencies(gamgee_test htslib)