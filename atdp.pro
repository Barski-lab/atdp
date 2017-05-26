######################################################################
# pro file for atdheatmap plot Wed Oct 19 14:59:16 2011
######################################################################

TEMPLATE = app
TARGET   = atdp

CONFIG   += console warn_on release
CONFIG   -= app_bundle

QT       -= gui

HEADERS     += src/Arguments.hpp \
               src/atdp.hpp \
               src/atdpbasics.hpp \
               src/Reads.hpp \
               src/main.hpp \
               src/config.hpp \
               src/math.hpp \
               src/Matrix.hpp \
               src/bam_reader_util.hpp

SOURCES     += src/Arguments.cpp \
               src/atdp.cpp \
               src/atdpbasics.cpp \
               src/Reads.cpp \
               src/main.cpp \
               src/bam_reader_util.cpp

INCLUDEPATH += . \
               ./src \
               ./bamtools

!win32{

OBJECTS_DIR = GeneratedFiles
UI_DIR      = GeneratedFiles
MOC_DIR     = GeneratedFiles
RCC_DIR     = GeneratedFiles

DEFINES        += _APPNAME=\\\"$$TARGET\\\"
LIBS           += -lm -lz ./bamtools/libbamtools.a -lz

lib_bamtools.commands = cd ./bamtools/; qmake; $(MAKE) -j 8
QMAKE_EXTRA_TARGETS   = lib_bamtools
PRE_TARGETDEPS        = lib_bamtools

INCLUDEPATH += /usr/local/include/
}

macx{
#QMAKE_CFLAGS_X86_64 += -mmacosx-version-min=10.7
QMAKE_CXXFLAGS_X86_64 = $$QMAKE_CFLAGS_X86_64
}

win32{
DEFINES        += _APPNAME=\"$$TARGET\"
LIBS           += -lbamtools
}

QMAKE_CLEAN += $${TARGET} logfile.log *~ ./src/*~ *.txt
