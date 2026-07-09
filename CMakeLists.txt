cmake_minimum_required(VERSION 3.15)
project(BHJet LANGUAGES CXX)

# --- Options ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# I added this only for unsuccesful profiling attempts
# set(CMAKE_BUILD_TYPE Debug)
# add_compile_options(-g)

# --- Dependencies ---
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(nanobind CONFIG REQUIRED)
find_package(GSL REQUIRED)

# Fetch Kariba from github
# include(FetchContent)
# FetchContent_Declare(
#   kariba
#   GIT_REPOSITORY https://github.com/maklinger/kariba.git
#   GIT_TAG bhjetnew-next
#   GIT_REPOSITORY https://github.com/antonpannekoek/kariba.git
#   GIT_REPOSITORY https://github.com/rubyduncan/kariba.git
#   GIT_TAG cutoff_new
# )
# FetchContent_MakeAvailable(kariba)

# # Collect Kariba sources from github
# file(GLOB KARIBA_SOURCES "${kariba_SOURCE_DIR}/src/*.cpp") #this is wrong? the wrong path
# include_directories(${kariba_SOURCE_DIR}/src)

set(kariba_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../kariba" CACHE PATH "Path to local Kariba source")
file(GLOB KARIBA_SOURCES "${kariba_SOURCE_DIR}/src/*.cpp")

add_library(kariba_lib STATIC ${KARIBA_SOURCES})
target_include_directories(kariba_lib PUBLIC
    ${kariba_SOURCE_DIR}/src
    ${GSL_INCLUDE_DIR}
)

# --- Sources ---
set(SRC_FILES
    src/BHJet.cpp
    src/BLJet.cpp
    src/IsoJet.cpp
    src/JetDynamics.cpp
    src/RadiationZone.cpp
    src/TargetBlackBody.cpp
    src/TargetPhotonField.cpp
    src/utils.cpp
    src/MLJet.cpp
)

# --- Python module ---
nanobind_add_module(bhjet
    litgen_bindings/module.cpp
    litgen_bindings/generated_bindings.cpp
    ${SRC_FILES}
)

target_compile_options(bhjet PRIVATE -Wall -Wextra -Wpedantic)


# Include bhjet headers
target_include_directories(bhjet PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link Kariba and other libraries
target_link_libraries(bhjet PRIVATE kariba_lib GSL::gsl GSL::gslcblas)
if(UNIX)
    target_link_libraries(bhjet PRIVATE m)
endif()


# --- Output ---
# installs library to /python/bhjet
set_target_properties(bhjet PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/bhjet
)


install(TARGETS bhjet DESTINATION python)


# --- Test executable ---
# add_executable(test_mljet
#     test/test.cpp
#     ${SRC_FILES}
# )

# target_include_directories(test_mljet PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}/include
# )

# target_link_libraries(test_mljet PRIVATE
#     kariba_lib
#     GSL::gsl
#     GSL::gslcblas
# )

# target_compile_options(test_mljet PRIVATE -Wall -Wextra -Wpedantic)

# if(UNIX)
#     target_link_libraries(test_mljet PRIVATE m)
# endif()
