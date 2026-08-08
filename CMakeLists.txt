cmake_minimum_required(VERSION 3.15)

project(ThermoRNA
    VERSION 0.1.0
    LANGUAGES CXX
    DESCRIPTION "ThermoRNA: A C++ library for processing thermodynamic data of RNA structures and their energy calculations."
)

# ------------------------- C++ Standard -------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------- Components -------------------------

option(THERMORNA_KNOTERGY "Build the Knotergy components" ON)
option(THERMORNA_VIENNARNA "Build the ViennaRNA components" ON)

if(NOT THERMORNA_KNOTERGY AND NOT THERMORNA_VIENNARNA)
    message(FATAL_ERROR "At least one of THERMORNA_KNOTERGY or THERMORNA_VIENNARNA must be enabled.")
endif()

message(STATUS "ThermoRNA components:")
message(STATUS "  Knotergy:  ${THERMORNA_KNOTERGY}")
message(STATUS "  ViennaRNA: ${THERMORNA_VIENNARNA}")

# ------------------------- Configuration -------------------------

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Enable ThermoRNA warnings when built directly.
# Disable them by default when ThermoRNA is consumed through FetchContent/add_subdirectory.
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(THERMORNA_WARNINGS_DEFAULT ON)
else()
    set(THERMORNA_WARNINGS_DEFAULT OFF)
endif()

option(THERMORNA_WARNINGS "Enable ThermoRNA compiler warnings" ${THERMORNA_WARNINGS_DEFAULT})

message(STATUS "  Warnings:   ${THERMORNA_WARNINGS}")

# ------------------------- External Dependencies -------------------------

include(FetchContent)

# nlohmann/json is only needed by the Knotergy subset.
if(THERMORNA_KNOTERGY)
    message(STATUS "Fetching external dependencies... this may take a while on the first run")

    set(JSON_Install OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.12.0
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# ------------------------- Source Files -------------------------

set(ThermoRNASources)

if(THERMORNA_KNOTERGY)
    file(GLOB_RECURSE ThermoRNAKnotergySources CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/Knotergy/*.cpp"
    )
    list(APPEND ThermoRNASources ${ThermoRNAKnotergySources})
endif()

if(THERMORNA_VIENNARNA)
    file(GLOB_RECURSE ThermoRNAViennaRNASources CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/ViennaRNA/*.cpp"
    )
    list(APPEND ThermoRNASources ${ThermoRNAViennaRNASources})
endif()

# ------------------------- ThermoRNA Library -------------------------

add_library(ThermoRNALib ${ThermoRNASources})
add_library(ThermoRNA::ThermoRNA ALIAS ThermoRNALib)

target_include_directories(ThermoRNALib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_compile_definitions(ThermoRNALib PRIVATE
    THERMORNA_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

if(THERMORNA_KNOTERGY)
    target_compile_definitions(ThermoRNALib PUBLIC THERMORNA_KNOTERGY=1)
    target_link_libraries(ThermoRNALib PUBLIC nlohmann_json::nlohmann_json)
endif()

if(THERMORNA_VIENNARNA)
    target_compile_definitions(ThermoRNALib PUBLIC THERMORNA_VIENNARNA=1)
endif()

# ------------------------- Compiler Flags -------------------------

if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "Strict")
    message(STATUS "Building ThermoRNA in Debug mode")

    target_compile_options(ThermoRNALib PRIVATE
        -g
        -O0
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
        -fno-omit-frame-pointer
    )

    if(THERMORNA_WARNINGS)
        target_compile_options(ThermoRNALib PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wshadow
            -Wsign-compare
            -Wconversion
            -Wsign-conversion
        )

        if(CMAKE_BUILD_TYPE STREQUAL "Strict")
            target_compile_options(ThermoRNALib PRIVATE -Werror)
        endif()
    endif()

    target_link_options(ThermoRNALib PRIVATE
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
    )
else()
    message(STATUS "Building ThermoRNA in Release/Normal mode")

    target_compile_options(ThermoRNALib PRIVATE -O3)

    if(THERMORNA_WARNINGS)
        target_compile_options(ThermoRNALib PRIVATE
            -Wall
            -Wextra
            -Wpedantic
        )
    endif()
endif()

# ------------------------- Install -------------------------

# install(TARGETS ThermoRNALib
#     ARCHIVE DESTINATION lib
#     LIBRARY DESTINATION lib
# )
#
# install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/"
#     DESTINATION include/ThermoRNA
#     FILES_MATCHING
#         PATTERN "*.hpp"
# )

# ------------------------- clang-format -------------------------

# Only create development/formatting targets when ThermoRNA is the top-level project.
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    find_program(CLANG_FORMAT_EXE NAMES clang-format)

    file(GLOB_RECURSE THERMORNA_FORMAT_FILES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/test/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/test/*.hpp"
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
endif()