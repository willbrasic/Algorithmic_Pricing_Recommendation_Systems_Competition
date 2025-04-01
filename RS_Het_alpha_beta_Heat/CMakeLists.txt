project(RS_Het_alpha_beta_Heat)

# Add the executable target
add_executable(RS_Het_alpha_beta_Heat main.cpp Functions.cpp Global.cpp)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Get compiler to optimize the most
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")