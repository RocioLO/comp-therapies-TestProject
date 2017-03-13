cmake_minimum_required(VERSION 3.6)

project(ctlab)

set(MITK_DIR "/home/lab/opt/apps/MITK-2016.03.0-src/build-release/MITK-build")
set(CMAKE_AUTOMOC ON)

find_package(MITK REQUIRED)
find_package(Qt5Widgets)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


mitk_create_executable(${PROJECT_NAME} DEPENDS MitkQtWidgetsExt)

