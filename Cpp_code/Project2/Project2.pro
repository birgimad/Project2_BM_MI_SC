QT += core
QT -= gui

TARGET = Project2
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp

LIBS += -llapack -lblas -larmadillo
