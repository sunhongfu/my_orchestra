cmake_minimum_required(VERSION 3.6)

project(my_orchestra)

# Uncomment for optional debugging
#set(CMAKE_VERBOSE_MAKEFILE ON)

# TOPDIR is the root of the source tree
<<<<<<< .merge_file_A168rl
set(OX_INSTALL_DIRECTORY "/Users/uqhsun8/code/orchestra-sdk-1.7-1")
=======
set(OX_INSTALL_DIRECTORY "/Users/hongfusun/code/orchestra-sdk-1.7-1")
>>>>>>> .merge_file_hQukmR
set(TOPDIR "${OX_INSTALL_DIRECTORY}/include")
set(LARGE_LIBRARY_TYPE STATIC)

message("Configuring GE SDK Recon Build")
set(PRODUCT EXTERNAL_SDK)

# Include recon libraries exported from the SDK CMake build
set(RECON_LIBRARIES_INCLUDE_FILE ${OX_INSTALL_DIRECTORY}/lib/ReconLibraries.cmake)
if(EXISTS ${RECON_LIBRARIES_INCLUDE_FILE})
    include (${RECON_LIBRARIES_INCLUDE_FILE})
else()
    message("Could not find ${RECON_LIBRARIES_INCLUDE_FILE}")
    message(FATAL_ERROR "Verify that the CMake OX_INSTALL_DIRECTORY option is set correctly")
endif()

# Include SDK build configuration
include (${TOPDIR}/recon/SDK/product.cmake)

# Include CMakeLists.txt for each rehearsal project
add_subdirectory (Calibration)
add_subdirectory (Cartesian)
add_subdirectory (Epi)
add_subdirectory (Flex)
add_subdirectory (Spectroscopy)
add_subdirectory (Spiral)

# downloaded projects
add_subdirectory (SWI)

# Hongfu's projects
add_subdirectory (HS_ModHeader)
add_subdirectory (HS_PrintHeader)