TEMPLATE = lib
CONFIG  += warn_off staticlib create_prl
#TARGET  += GeometryFileTypes

macx-g++ {
  QMAKE_CFLAGS += -m64
  QMAKE_CXXFLAGS += -m64
}


# Input
INCLUDEPATH += ../Geometry

SOURCES =  \
		GeometryFileType.cpp \
		GeometryLoader.cpp \
		MayaOBJFile.cpp \
		RawcFile.cpp \
		RawFile.cpp \
		RawncFile.cpp \
		RawnFile.cpp


HEADERS =  \
		GeometryFileType.h \
		GeometryLoader.h \
		MayaOBJFile.h \
		RawcFile.h \
		RawFile.h \
		RawncFile.h \
		RawnFile.h


