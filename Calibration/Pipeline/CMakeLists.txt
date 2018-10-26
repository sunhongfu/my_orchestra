project(CalibrationRehearsalPipeline)

set(SOURCE_FILES
        Calibration2DRecon.cpp
        Calibration3DRecon.cpp
        CalibrationRehearsalUtils.cpp
        Calibration2DRecon.h
        Calibration3DRecon.h
        CalibrationRehearsalUtils.h)

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES})

# Install this example rehearsal code along with this CMakeLists.txt file
install(FILES ${SOURCE_FILES} DESTINATION "Examples/Calibration/Pipeline")
install(FILES "CMakeLists.txt" DESTINATION "Examples/Calibration/Pipeline")
