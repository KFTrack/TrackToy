add_library(Detector SHARED
    Calorimeter.cc
    Tracker.cc
    Target.cc
    HollowCylinder.cc
    IPA.cc
    EStar.cc
    TrckrWall.cc
    )

# set top-level directory as include root
target_include_directories(Detector PRIVATE ${PROJECT_SOURCE_DIR}/..)
# add KinKal; this may be redundant, but...
target_include_directories(Detector PRIVATE ${KINKAL_INCLUDE_DIR})
# add root
target_link_libraries(Detector ${ROOT_LIBRARIES} )

# set shared library version equal to project version
set_target_properties(Detector PROPERTIES VERSION ${PROJECT_VERSION})
