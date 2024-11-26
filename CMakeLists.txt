cmake_minimum_required(VERSION 3.15...3.27)

project(dehasher LANGUAGES CXX)

if (NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()

# Find Python components
find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

# Find nanobind
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Add subdirectories for additional libraries if needed
add_subdirectory(lib/json)
add_subdirectory(lib/kseqpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/simdjson)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/json/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/includes)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/kseqpp/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/parallel-hashmap)


# Define the nanobind module
nanobind_add_module(
  _dehasher_impl
  STABLE_ABI
  NB_STATIC
  src/hashutil.cpp
  src/dehasher.cpp 
  lib/simdjson/simdjson.cpp
)

# Link libraries properly
target_link_libraries(_dehasher_impl PRIVATE 
  nlohmann_json::nlohmann_json 
  kseq++::kseq++
)

# Installation
install(TARGETS _dehasher_impl LIBRARY DESTINATION dehasher)