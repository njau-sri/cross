TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_LFLAGS += -static

SOURCES += main.cpp \
    appcross.cpp \
    cmdline.cpp \
    util.cpp \
    hapmapio.cpp \
    plinkio.cpp \
    rtmio.cpp \
    vcfio.cpp

HEADERS += \
    appcross.h \
    main.h \
    cmdline.h \
    util.h \
    strsplit.h \
    hapmapio.h \
    plinkio.h \
    rtmio.h \
    vcfio.h
