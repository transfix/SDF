TEMPLATE = lib
CONFIG  += warn_off staticlib create_prl
#TARGET  += GeometryFileTypes

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


