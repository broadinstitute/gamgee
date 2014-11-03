# build htslib
set(htslib_PREFIX ${CMAKE_BINARY_DIR}/contrib/htslib)
ExternalProject_Add(htslib
    PREFIX ${htslib_PREFIX}
    GIT_REPOSITORY "https://github.com/broadinstitute/htslib.git"
    GIT_TAG broad
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make lib-static -j 4
    INSTALL_COMMAND ""
    LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
)

include_directories(${htslib_PREFIX}/src/htslib)
set(htslib_LIB ${htslib_PREFIX}/src/htslib/libhts.a)