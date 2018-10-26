add_subdirectory(Pipeline)

project(CalibrationRehearsal)

set(SOURCE_FILES
        Rehearsal.cpp
        Rehearsal.h)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} CalibrationRehearsalPipeline)
target_link_libraries(${PROJECT_NAME} Calibration3D)
target_link_libraries(${PROJECT_NAME} CalibrationCommon)
target_link_libraries(${PROJECT_NAME} Pure1)
target_link_libraries(${PROJECT_NAME} Clariview)
target_link_libraries(${PROJECT_NAME} Asset)
target_link_libraries(${PROJECT_NAME} SpectroMCSI)
target_link_libraries(${PROJECT_NAME} SpectroMCSILegacy)
target_link_libraries(${PROJECT_NAME} Cartesian3D)
target_link_libraries(${PROJECT_NAME} Cartesian2D)
target_link_libraries(${PROJECT_NAME} Gradwarp)
target_link_libraries(${PROJECT_NAME} Legacy)
target_link_libraries(${PROJECT_NAME} Core)
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
install(FILES ${SOURCE_FILES} DESTINATION "Examples/Calibration")
install(FILES "CMakeLists.txt" DESTINATION "Examples/Calibration")
