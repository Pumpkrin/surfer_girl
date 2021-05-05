cmake_minimum_required(VERSION 3.2)

project(surfer_girl_project VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT 6.10 CONFIG REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(RootNewMacros)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib/)

root_generate_dictionary( D__surfer_girl data_structure.hpp LINKDEF linkdef.h)

add_library( data_structure SHARED D__surfer_girl.cxx )
target_include_directories( data_structure PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries( data_structure PUBLIC ROOT::Hist)
set_target_properties( data_structure PROPERTIES SUFFIX ".so" )

add_executable( surfer_girl surfer_girl.cpp )
target_include_directories(surfer_girl PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(surfer_girl PUBLIC ROOT::Core ROOT::Tree ROOT::Hist data_structure)

add_executable( tree_flip tree_flip.cpp )
target_include_directories(tree_flip PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(tree_flip PUBLIC ROOT::Core ROOT::Tree ROOT::Hist data_structure)

add_executable( tree_grubbing tree_grubbing.cpp )
target_include_directories(tree_grubbing PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(tree_grubbing PUBLIC ROOT::Core ROOT::Tree ROOT::Hist data_structure)

set( path_to_lib "${CMAKE_INSTALL_LIBDIR}") 
set( path_to_macro "${CMAKE_BINARY_DIR}/macros/") 
configure_file( rootlogon.C rootlogon.C @ONLY) 
configure_file( reprocess_waveform.hpp macros/reprocess_waveform.hpp COPYONLY ) 

