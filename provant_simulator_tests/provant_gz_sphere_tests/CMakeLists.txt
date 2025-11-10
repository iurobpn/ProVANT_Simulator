# This file is part of the ProVANT simulator project. Licensed under the terms of the MIT open source license. More
# details at https: //github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md

cmake_minimum_required(VERSION 3.0.2)
project(provant_gz_sphere_tests)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.16.0 # Or another desired version/tag
    )
FetchContent_MakeAvailable(spdlog)
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 11.2.0 # Or a more recent version tag
)
FetchContent_MakeAvailable(fmt)
find_package(
  catkin REQUIRED
  COMPONENTS gazebo_ros
             provant_simulator_cmake
             provant_simulator_log_utils
             provant_simulator_sdf_parser
             roscpp
)
find_package(gazebo REQUIRED)
find_package(Boost REQUIRED)

set(PLUGIN_LIBRARIES
    provant_test_apply_force_plugin
    provant_test_apply_rel_force_once_plugin
    provant_test_set_force_once_plugin
    provant_test_set_force_plugin
)

catkin_package(
  INCLUDE_DIRS
  include
  LIBRARIES
  ${PLUGIN_LIBRARIES}
  CATKIN_DEPENDS
  gazebo_ros
  provant_simulator_cmake
  provant_simulator_log_utils
  provant_simulator_sdf_parser
  roscpp
  DEPENDS
  Boost
)

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
  ${GAZEBO_INCLUDE_DIRS}
  ${fmt_INCLUDE_DIRS}
  include
)

add_library(provant_test_apply_rel_force_once_plugin src/apply_force_once/apply_force_once.cpp)
add_dependencies(
  provant_test_apply_rel_force_once_plugin ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS}
)
target_link_libraries(
  provant_test_apply_rel_force_once_plugin
  ${catkin_LIBRARIES}
  ${GAZEBO_LIBRARIES}
  ${Boost_LIBRARIES}
  spdlog::spdlog
)

add_library(provant_test_set_force_once_plugin src/set_force_once/set_force_once.cpp)
add_dependencies(provant_test_set_force_once_plugin ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(
  provant_test_set_force_once_plugin
  ${catkin_LIBRARIES}
  ${GAZEBO_LIBRARIES}
  ${Boost_LIBRARIES}
  spdlog::spdlog
)

add_library(provant_test_apply_force_plugin src/apply_force_multiple/apply_force_multiple.cpp)
add_dependencies(provant_test_apply_force_plugin ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(
  provant_test_apply_force_plugin
  ${catkin_LIBRARIES}
  ${GAZEBO_LIBRARIES}
  ${Boost_LIBRARIES}
  spdlog::spdlog
)

add_library(provant_test_set_force_plugin src/set_force/set_force.cpp)
add_dependencies(provant_test_set_force_plugin ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(
  provant_test_set_force_plugin
  ${catkin_LIBRARIES}
  ${GAZEBO_LIBRARIES}
  ${Boost_LIBRARIES}
  spdlog::spdlog
)

install(
  TARGETS ${PLUGIN_LIBRARIES}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)

install(
  DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN ".svn" EXCLUDE
)
