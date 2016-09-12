# For private use only.
# Copyright (c) Manuel Freiberger, 2016.

TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

HEADERS += \
    ../GaloisField.hpp \
    ../ReedSolomon.hpp

SOURCES += main.cpp \
    tst_gf256element.cpp \
    tst_reedsolomon.cpp

