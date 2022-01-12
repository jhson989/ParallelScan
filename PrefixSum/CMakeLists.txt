# CMake bianry version
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

set(APP_NAME app.out)
set(MAIN_FILE main.cpp)
set(COMPILE_OPTIONS -fopenmp)

add_executable (${APP_NAME} ${MAIN_FILE})
target_link_options(${APP_NAME} PRIVATE ${COMPILE_OPTIONS})