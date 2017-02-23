cmake_minimum_required(VERSION 3.6)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(test)

INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(/usr/include/eigen3)

FILE(
    GLOB_RECURSE
    source_files
    include/*
    src/*
)

#ADD_EXECUTABLE(main test/main.cpp ${source_files})
ADD_EXECUTABLE(main_cartpole test/main_cartpole.cpp ${source_files})
#ADD_EXECUTABLE(mainMPC test/mainMPC.cpp ${source_files})
#ADD_EXECUTABLE(testModel test/testModel.cpp ${source_files})

#TARGET_LINK_LIBRARIES(main ${qpOASES_LIBRARIES})
TARGET_LINK_LIBRARIES(main_cartpole ${qpOASES_LIBRARIES})
#TARGET_LINK_LIBRARIES(mainMPC ${qpOASES_LIBRARIES})
#TARGET_LINK_LIBRARIES(testModel ${qpOASES_LIBRARIES})
