cmake_minimum_required(VERSION 3.14)

# ---- Project ----

# Note: update this to your new project's name and version
project(
  kseq
  VERSION 1.31
  LANGUAGES CXX
)

# ---- Include CPM ----
include(${PROJECT_SOURCE_DIR}/cmake/CPM.cmake)


# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()


# PackageProject.cmake will be used to make our target installable
CPMAddPackage(
  NAME PackageProject.cmake
  GITHUB_REPOSITORY TheLartians/PackageProject.cmake
  VERSION 1.5.0
)

file(GLOB_RECURSE headers CONFIGURE_DEPENDS include/*.h)

# ---- Create library ----
add_library(kseq INTERFACE)

target_include_directories(
  kseq INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                 $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)


string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)


packageProject(
  NAME ${PROJECT_NAME}
  VERSION ${PROJECT_VERSION}
  NAMESPACE ${PROJECT_NAME}
  BINARY_DIR ${PROJECT_BINARY_DIR}
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
  INCLUDE_DESTINATION include/${PROJECT_NAME}-${PROJECT_VERSION}
  VERSION_HEADER "${VERSION_HEADER_LOCATION}"
  COMPATIBILITY SameMajorVersion
  DEPENDENCIES ""
)