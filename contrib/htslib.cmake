if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # when using the makefile generator, use the special variable $(MAKE) to invoke make
    # this enables the jobserver to work correctly
    set(MAKE_COMMAND "$(MAKE)")
else()
    # invoke make explicitly
    # in this case, we assume the parent build system is running in parallel already so no -j flag is added
    find_program(MAKE_COMMAND NAMES make gmake)
endif()

if (INSTALL_DEPENDENCIES)
    set(HTSLIB_INSTALL ${MAKE_COMMAND} install prefix=${CMAKE_INSTALL_PREFIX})
else()
    set(HTSLIB_INSTALL "")
endif()

# build htslib
set(htslib_PREFIX ${CMAKE_BINARY_DIR}/contrib/htslib)
ExternalProject_Add(htslib
    PREFIX ${htslib_PREFIX}
    GIT_REPOSITORY "https://github.com/samtools/htslib.git"
    GIT_TAG 49fdfbda20acbd73303df3c7fef84f2d972c5f8d # Version 1.5
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${MAKE_COMMAND} lib-static
    INSTALL_COMMAND "${HTSLIB_INSTALL}"
    LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
)

include_directories(${htslib_PREFIX}/src/htslib)
set(htslib_LIB ${htslib_PREFIX}/src/htslib/libhts.a)
