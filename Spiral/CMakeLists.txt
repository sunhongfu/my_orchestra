project(SpiralRehearsal)

set(SOURCE_FILES
    Rehearsal.cpp
    Spiral2D.cpp
    SpiralScanArchiveRecon.cpp
    SpiralUtil.cpp
    Rehearsal.h
    SpiralUtil.h)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} Arc)
target_link_libraries(${PROJECT_NAME} Cartesian2D)
target_link_libraries(${PROJECT_NAME} Cartesian3D)
target_link_libraries(${PROJECT_NAME} Gradwarp)
target_link_libraries(${PROJECT_NAME} Spiral)
target_link_libraries(${PROJECT_NAME} Legacy)
target_link_libraries(${PROJECT_NAME} Core)
target_link_libraries(${PROJECT_NAME} CalibrationCommon)
target_link_libraries(${PROJECT_NAME} Foundation)
target_link_libraries(${PROJECT_NAME} Acquisition)
target_link_libraries(${PROJECT_NAME} Control)
target_link_libraries(${PROJECT_NAME} Common)
target_link_libraries(${PROJECT_NAME} Crucial)
target_link_libraries(${PROJECT_NAME} Dicom)
target_link_libraries(${PROJECT_NAME} Hdf5)
target_link_libraries(${PROJECT_NAME} Math)
target_link_libraries(${PROJECT_NAME} SystemServicesImplementation)
target_link_libraries(${PROJECT_NAME} SystemServicesInterface)
target_link_libraries(${PROJECT_NAME} System)
target_link_libraries(${PROJECT_NAME} ${OX_3P_LIBS})
target_link_libraries(${PROJECT_NAME} ${OX_OS_LIBS})

# Install this example rehearsal code along with this CMakeLists.txt file
install(FILES ${SOURCE_FILES} DESTINATION "Examples/Spiral")
install(FILES "CMakeLists.txt" DESTINATION "Examples/Spiral")
