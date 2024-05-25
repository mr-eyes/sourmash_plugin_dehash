# kseq (CMake Friendly)

[![Test_Include](https://github.com/mr-eyes/kseq/actions/workflows/standalone.yml/badge.svg?branch=master)](https://github.com/mr-eyes/kseq/actions/workflows/standalone.yml)

This is just a CMake wrapper for [kseq.h](http://lh3lh3.users.sourceforge.net/kseq.shtml).

## Usage

### Include

```cpp
#include <kseq/kseq.h>
```

#### Retreive the whole header comment

In Fasta headers, `kseq.h` can return two parts.

1. `kseqObj->name.s`: The sequence id, Or technincally speaking, the first word before the delimiter (`\t`, `\s`).
2. `kseqObj->comment.s`: The rest of the header after the delimiter.

To retreive the whole header you would do something like: `kseqObj->name.s + delimiter + kseqObj->comment.s`. The delimiter here is unknown due to its variablity in the sequence files.

To get the full comment, just put the following line inside the main function, or your target function.

```cpp
KS_FULL_COMMENT = true; 
```

Now the `kseqObj->comment.s` = `delimiter` + the rest of the header line.

### Method 1 (add_subdirectory)

```cmake
add_subdirectory(kseq)
add_executable(ex1 main.cpp)
target_link_libraries(ex1 kseq)
```

### Method 2 ([CPM.cmake](https://github.com/cpm-cmake/CPM.cmake))

```cmake
# ---- Setting up CPM.cmake

include(${PROJECT_SOURCE_DIR}/CPM.cmake)

# ---- Including kseq library

## Online?
CPMAddPackage(
  NAME kseq
  GIT_TAG 1.3
  GITHUB_REPOSITORY mr-eyes/kseq # to get an installable target)
)

## Have the source Offline?
if(FALSE) # fake a block comment
CPMAddPackage(
  NAME kseq
  SOURCE_DIR  ${CMAKE_CURRENT_LIST_DIR}/kseq # to get an installable target)
)
endif()


add_executable(ex2 main.cpp)
target_link_libraries(ex2 kseq)
```
