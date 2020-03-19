SOURCES += main.cpp \
    main_backup.cpp
HEADERS += src/*.h

INCLUDEPATH += \
            $$PWD/src \
            /home/shifeng/Codes/0.Codes \
            /usr/lib/gcc

TARGET = main
LIBS += -L/opt/intel/mkl/lib/intel64
