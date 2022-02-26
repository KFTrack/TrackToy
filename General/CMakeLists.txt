
add_library(General SHARED
    MuonRange.cc
    ELossDistributions.cc
)

# set top-level directory as include root
target_include_directories(General PRIVATE ${PROJECT_SOURCE_DIR}/.. ${KINKAL_INCLUDE_DIR})

# set shared library version equal to project version
set_target_properties(General PROPERTIES VERSION ${PROJECT_VERSION})
