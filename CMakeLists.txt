cmake_minimum_required(VERSION 3.15)

project(ThermoRNA
    LANGUAGES CXX
    DESCRIPTION "ThermoRNA: A C++ library for processing thermodynamic data of RNA structures and their energy calculations."
)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Define a global compile definition for the source directory,
# so it can be used in both the library and tests
add_compile_definitions(
    THERMORNA_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

# ------------------------- External Dependencies ------------------------- 
include(FetchContent)

message(STATUS "Fetching external dependencies... this may take a while on the first run")

# ----------------- nlohmann/json -----------------
set(JSON_Install OFF CACHE BOOL "" FORCE)

FetchContent_Declare(nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.12.0
)
FetchContent_MakeAvailable(nlohmann_json)

# ------------------------- Configurations -------------------------

# Default build type if none is specified

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# -------------------------------------------------------------------
#                             C++ Flags
# -------------------------------------------------------------------

# Debug configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "Strict")
    message(STATUS "Building in Debug mode")
    add_compile_options(
        -g
        -O0
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wsign-compare
        -Wconversion
        -Wsign-conversion
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
        -fno-omit-frame-pointer
    )

    if(CMAKE_BUILD_TYPE STREQUAL "Strict")
        add_compile_options(-Werror)
    endif()

    add_link_options(
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
    )
else()
    # Normal (Release or other non-Debug builds)
    message(STATUS "Building in Release/Normal mode")
    add_compile_options(
        -O3
        -Wall
        -Wextra
        -Wpedantic
    )
endif()


# # -------------------------------------------------------------------
# #                      GoogleTest Configuration
# # -------------------------------------------------------------------
# option(BUILD_TESTING "Build unit tests" ON)

# if (NOT BUILD_TESTING)
#     message(STATUS "Unit tests are disabled. To enable them, set -DBUILD_TESTING=ON")
# elseif(EXISTS "${CMAKE_SOURCE_DIR}/test/googletest/CMakeLists.txt")
#     enable_testing()
#     add_subdirectory(test)
#     add_subdirectory(test/googletest)
# else()
#     message(WARNING "
# ============================================================
#  GoogleTest submodule not found. Unit tests will be disabled.
# ------------------------------------------------------------
# To enable unit testing, please do one of the following:

#   1. If you already cloned the repository:
#      Run the following command to initialize submodules:
#          git submodule update --init --recursive

#   2. If you haven't cloned the repository yet:
#      Clone it with submodules using:
#          git clone --recursive <repo-url>

# After initializing the submodule, re-run CMake to enable tests.
# ============================================================")
# endif()


# -------------------------------------------------------------------
#                     ThermoRNA Configuration
# -------------------------------------------------------------------

# ------------------------- Source Files -------------------------
file(GLOB_RECURSE ThermoRNASources
    ${CMAKE_SOURCE_DIR}/src/*.cpp
)

add_library(ThermoRNALib ${ThermoRNASources})

target_include_directories(ThermoRNALib PUBLIC ${CMAKE_SOURCE_DIR}/src)


# ------------------------- Link Libraries -------------------------
target_link_libraries(ThermoRNALib
  PUBLIC
    nlohmann_json::nlohmann_json
)

# # ------------------------- Install -------------------------
# install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
# install(DIRECTORY ${CMAKE_SOURCE_DIR}/params DESTINATION share/thermorna)

# -------------------------------------------------------------------
#                         clang-format
# -------------------------------------------------------------------

find_program(CLANG_FORMAT_EXE NAMES clang-format)

file(GLOB_RECURSE THERMORNA_FORMAT_FILES
    CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_SOURCE_DIR}/test/*.cpp"
    "${CMAKE_SOURCE_DIR}/test/*.hpp"
)

if(CLANG_FORMAT_EXE)
    add_custom_target(format
        COMMAND ${CLANG_FORMAT_EXE} -i ${THERMORNA_FORMAT_FILES}
        COMMENT "Formatting C++ source files with clang-format"
        VERBATIM
    )

    add_custom_target(format-check
        COMMAND ${CLANG_FORMAT_EXE} --dry-run --Werror ${THERMORNA_FORMAT_FILES}
        COMMENT "Checking C++ source formatting with clang-format"
        VERBATIM
    )
else()
    add_custom_target(format
        COMMAND ${CMAKE_COMMAND} -E echo "clang-format not found"
        COMMAND ${CMAKE_COMMAND} -E false
        VERBATIM
    )

    add_custom_target(format-check
        COMMAND ${CMAKE_COMMAND} -E echo "clang-format not found"
        COMMAND ${CMAKE_COMMAND} -E false
        VERBATIM
    )
endif()