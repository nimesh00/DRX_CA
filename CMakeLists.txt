cmake_minimum_required(VERSION 3.10)
project(DRX_CA CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to a Release build if the user didn't specify one.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# SDL2 — prefer the upstream CMake config, fall back to pkg-config.
find_package(SDL2 QUIET CONFIG)
if(NOT SDL2_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(SDL2 REQUIRED IMPORTED_TARGET sdl2)
    set(SDL2_TARGET PkgConfig::SDL2)
else()
    if(TARGET SDL2::SDL2)
        set(SDL2_TARGET SDL2::SDL2)
    else()
        # Older SDL2 CMake configs only export include dirs + libs as variables.
        set(SDL2_TARGET "")
    endif()
endif()

add_executable(drx
    drx.cpp
    viewer.cpp
)

# viewer.cpp is the only TU that needs SDL2; drx.cpp only sees viewer.h.
if(SDL2_TARGET)
    target_link_libraries(drx PRIVATE ${SDL2_TARGET})
else()
    target_include_directories(drx PRIVATE ${SDL2_INCLUDE_DIRS})
    target_link_libraries(drx PRIVATE ${SDL2_LIBRARIES})
endif()

# Threading for std::this_thread::sleep_for.
find_package(Threads REQUIRED)
target_link_libraries(drx PRIVATE Threads::Threads)
