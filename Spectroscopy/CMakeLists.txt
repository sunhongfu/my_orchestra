project(SpectroRehearsal)

set(SOURCE_FILES
    Rehearsal.cpp
    MultiChannelMultiVoxel.cpp
    PerChannelMultiVoxel.cpp
    SingleVoxel.cpp
    MultiChannelMultiVoxel.h
    PerChannelMultiVoxel.h
    Rehearsal.h
    SingleVoxel.h)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} SpectroMCSI)
target_link_libraries(${PROJECT_NAME} SpectroMCSILegacy)
target_link_libraries(${PROJECT_NAME} SpectroMultiVoxel)
target_link_libraries(${PROJECT_NAME} SpectroSingleVoxel)
target_link_libraries(${PROJECT_NAME} SpectroCommon)
target_link_libraries(${PROJECT_NAME} Legacy)
target_link_libraries(${PROJECT_NAME} Core)
target_link_libraries(${PROJECT_NAME} CalibrationCommon)
target_link_libraries(${PROJECT_NAME} Control)
target_link_libraries(${PROJECT_NAME} Common)
target_link_libraries(${PROJECT_NAME} Crucial)
target_link_libraries(${PROJECT_NAME} Dicom)
target_link_libraries(${PROJECT_NAME} Hdf5)
target_link_libraries(${PROJECT_NAME} Math)
target_link_libraries(${PROJECT_NAME} System)
target_link_libraries(${PROJECT_NAME} ${OX_3P_LIBS})
target_link_libraries(${PROJECT_NAME} ${OX_OS_LIBS})

# Install this example rehearsal code along with this CMakeLists.txt file
install(FILES ${SOURCE_FILES} DESTINATION "Examples/Spectroscopy")
install(FILES "CMakeLists.txt" DESTINATION "Examples/Spectroscopy")
