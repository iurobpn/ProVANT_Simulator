# This file is part of the ProVANT simulator project. Licensed under the terms of the MIT open source license. More
# details at https://github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md

cmake_minimum_required(VERSION 3.0.2)
project(vant4_2_winfnonlinear CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(cmake_modules REQUIRED)
find_package(
  catkin REQUIRED
  COMPONENTS roscpp
             std_msgs
             simulator_msgs
             control_strategies_base
             provant_simulator_cmake
)
find_package(
  Eigen3
  3.3
  REQUIRED
  NO_MODULE
)
find_package(CPLEX)

if(CPLEX_FOUND)
  catkin_package(
    INCLUDE_DIRS
    include
    LIBRARIES
    ${PROJECT_NAME}
    CATKIN_DEPENDS
    roscpp
    std_msgs
    simulator_msgs
    control_strategies_base
    # DEPENDS
  )

  include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${CPLEX_INCLUDE_DIR}
    ${CPLEX_CONCERT_INCLUDE_DIR}
  )

  add_library(${PROJECT_NAME} src/cplexutils.cpp src/main.cpp)
  add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
  target_link_libraries(
    ${PROJECT_NAME}
    ${catkin_LIBRARIES}
    ${CPLEX_LIBRARIES}
    ${CPLEX_CONCERT_LIBRARY}
    ${CPLEX_EXTRA_LIBRARIES}
    Eigen3::Eigen
  )
  add_compile_definitions(${CPLEX_COMPILER_DEFINITIONS})

  # Padrao do simulador
  install(
    TARGETS ${PROJECT_NAME}
    ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )
  install(
    DIRECTORY include
    DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
    PATTERN ".svn" EXCLUDE
  )

else()
  catkin_package()
  message(WARNING "CPLEX was not found, therefore the ${PROJECT_NAME} will not be compiled.")
endif()
