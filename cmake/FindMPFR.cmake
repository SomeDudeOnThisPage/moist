find_path(MPFR_INCLUDE_DIR
    NAMES mpfr.h
    PATHS /usr/include /usr/local/include /opt/homebrew/include
)

find_library(MPFR_LIBRARIES
    NAMES mpfr
    PATHS /usr/lib /usr/local/lib /opt/homebrew/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_LIBRARIES MPFR_INCLUDE_DIR)

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARIES)
