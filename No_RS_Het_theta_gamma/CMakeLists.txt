cmake_minimum_required(VERSION 3.26.5)

project(No_RS_Het_theta_gamma)

# Add the executable target
add_executable(No_RS_Het_theta_gamma main.cpp Functions.cpp Global.cpp)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Get compiler to optimize the most
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")