TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    cmdline.cpp \
    cross.cpp \
    vcfio.cpp

HEADERS += \
    cmdline.h \
    split.h \
    util.h \
    vcfio.h
