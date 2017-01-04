cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set(APP_NAME "Evolve")
project(${APP_NAME})

set(PROJECT_LIB_DIR  ${${APP_NAME}__SOURCE_DIR}../Common/lib)
set(PROJECT_INC_DIR  ${${APP_NAME}__SOURCE_DIR}../Common)

## OpenMesh
find_package (OpenMesh)
if (NOT OPENMESH_FOUND)
    message (FATAL_ERROR "OpenMesh not found!")
endif ()


## Eigen 
SET( EIGEN3_INCLUDE_DIR "D:/Program Files/eigen-eigen-04ab5fa4b241" )
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )



FILE(GLOB SRC_SOURCES *.cpp)
FILE(GLOB SRC_HEADERS *.h )


set(SRC_LIST ${SRC_SOURCES} 
              ${SRC_HEADERS}
              ${SRC_MATH}
              )
              
include_directories(
  ${PROJECT_INC_DIR}
  ${OPENMESH_INCLUDE_DIRS}
)
link_directories(
  ${PROJECT_LIB_DIR}
  ${OPENMESH_LIBRARIES}
)              
              
######################################################################			      
set(CMAKE_BUILD_TYPE Release)

add_executable (${APP_NAME} ${SRC_LIST} )

target_link_libraries (${APP_NAME}
    ${OPENMESH_LIBRARIES}    	
)
