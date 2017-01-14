TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    src/appcross.cpp \
    src/cmdline.cpp \
    src/hapmapio.cpp \
    src/main.cpp \
    src/plinkio.cpp \
    src/rtmio.cpp \
    src/util.cpp \
    src/vcfio.cpp

HEADERS += \
    src/appcross.h \
    src/cmdline.h \
    src/hapmapio.h \
    src/main.h \
    src/plinkio.h \
    src/rtmio.h \
    src/strsplit.h \
    src/util.h \
    src/vcfio.h
