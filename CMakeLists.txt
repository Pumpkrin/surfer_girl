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

add_library( fit_tools INTERFACE )
target_include_directories( fit_tools INTERFACE "${ROOT_INCLUDE_DIRS}")
add_library( formulae INTERFACE )
target_include_directories( formulae INTERFACE "${ROOT_INCLUDE_DIRS}")
target_link_libraries( formulae INTERFACE fit_tools)

add_executable( surfer_girl surfer_girl.cpp )
target_include_directories(surfer_girl PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(surfer_girl PUBLIC ROOT::Core ROOT::Tree ROOT::Hist data_structure)

add_executable( tree_flip tree_flip.cpp )
target_include_directories(tree_flip PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(tree_flip PUBLIC ROOT::Core ROOT::Tree ROOT::Hist data_structure)

add_executable( tree_grubbing tree_grubbing.cpp )
target_include_directories(tree_grubbing PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(tree_grubbing PUBLIC ROOT::Core ROOT::Tree ROOT::Hist data_structure)

add_executable( energy_loss energy_loss.cpp )
target_include_directories(energy_loss PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(energy_loss PUBLIC ROOT::Core ROOT::Graf)

add_executable( calibration  calibration.cpp )
target_include_directories(calibration PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(calibration PUBLIC ROOT::Core ROOT::Graf fit_tools )

add_executable( energy_resolution  energy_resolution.cpp )
target_include_directories(energy_resolution PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(energy_resolution PUBLIC ROOT::Core ROOT::Graf fit_tools )

add_executable( generate_response_matrix  generate_response_matrix.cpp )
target_include_directories(generate_response_matrix PUBLIC "${ROOT_INCLUDE_DIRS}")
target_link_libraries(generate_response_matrix PUBLIC ROOT::Core ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad data_structure )

set( path_to_lib "${CMAKE_INSTALL_LIBDIR}") 
set( path_to_macro "${CMAKE_BINARY_DIR}/macros/") 
configure_file( rootlogon.C rootlogon.C @ONLY) 
configure_file( reprocess_waveform.hpp macros/reprocess_waveform.hpp COPYONLY ) 
configure_file( tools.hpp macros/tools.hpp COPYONLY ) 
configure_file( fit_tools.hpp macros/fit_tools.hpp COPYONLY ) 
configure_file( formulae.hpp macros/formulae.hpp COPYONLY ) 
configure_file( process_spectra.hpp macros/process_spectra.hpp COPYONLY ) 
configure_file( process_directory.sh macros/process_directory.sh COPYONLY ) 
configure_file( process_file.sh macros/process_file.sh COPYONLY ) 
configure_file( apply_cut.sh macros/apply_cut.sh COPYONLY ) 
configure_file( data/p_aluminum.txt data/p_aluminum.txt COPYONLY ) 
configure_file( data/p_air.txt data/p_air.txt COPYONLY ) 
configure_file( data/p_glass.txt data/p_glass.txt COPYONLY ) 
configure_file( data/p_teflon.txt data/p_teflon.txt COPYONLY ) 

