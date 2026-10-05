cmake_minimum_required(VERSION 3.16)
project(quantum_stormsynth VERSION 0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(QSTORM_WARNINGS -Wall -Wextra -Wpedantic)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FLAC REQUIRED IMPORTED_TARGET flac)
find_package(GTest REQUIRED)
set(QSTORM_QUERY_RUNNER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/qr")
add_subdirectory(${QSTORM_QUERY_RUNNER_DIR} ${CMAKE_BINARY_DIR}/query_runner)
file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS
  src/atom/*.cpp
  src/io/*.cpp
  src/cli/*.cpp
  src/ui/*.cpp
  src/*.cpp
  )

list(FILTER CORE_SOURCES EXCLUDE REGEX "src/main\\.cpp$")
set(QSTORM_GENERATED_DIR ${CMAKE_BINARY_DIR}/generated)
add_custom_command(
  OUTPUT ${QSTORM_GENERATED_DIR}/atom_version.h
  COMMAND ${CMAKE_COMMAND} -DSOURCE=${CMAKE_SOURCE_DIR} -DOUTPUT=${QSTORM_GENERATED_DIR}/atom_version.h -P ${CMAKE_SOURCE_DIR}/cmake/atom_version.cmake
  DEPENDS ${CMAKE_SOURCE_DIR}/src/atom/atom_thunder.h ${CMAKE_SOURCE_DIR}/src/atom/atom_thunder.cpp ${CMAKE_SOURCE_DIR}/.git/HEAD ${CMAKE_SOURCE_DIR}/cmake/atom_version.cmake
)
list(APPEND CORE_SOURCES ${QSTORM_GENERATED_DIR}/atom_version.h)
add_library(qstorm_core STATIC ${CORE_SOURCES})
target_include_directories(qstorm_core PUBLIC ${CMAKE_SOURCE_DIR}/src ${QSTORM_GENERATED_DIR})
target_compile_options(qstorm_core PRIVATE ${QSTORM_WARNINGS})
target_link_libraries(qstorm_core PUBLIC
  PkgConfig::FLAC
  query_runner
)
target_compile_definitions(qstorm_core PUBLIC
  QSTORM_HAS_FLAC
  QSTORM_HAS_PQXX
)

add_executable(qstorm src/main.cpp)
target_link_libraries(qstorm PRIVATE qstorm_core)
target_compile_options(qstorm PRIVATE ${QSTORM_WARNINGS})

enable_testing()
include(GoogleTest)
