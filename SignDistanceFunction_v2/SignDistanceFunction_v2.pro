TEMPLATE = lib
CONFIG  += qt warn_off opengl staticlib create_prl
#TARGET  += SignDistanceFunction_v2

macx-g++ {
	INCLUDEPATH += /usr/X11R6/include
}

# Input

SOURCES =  \
		bufferedio.cpp \
		DistanceTransform.cpp \
		FaceVertSet3D.cpp \
		geom.cpp \
		Geom3DParser.cpp \
		mtxlib.cpp \
		RawivParser.cpp \
		Reg3Parser.cpp


HEADERS =  \
		bio.h \
		bufferedio.h \
		cubes.h \
		diskio.h \
		DistanceTransform.h \
		dynarray.h \
		FaceVertSet3D.h \
		geom.h \
		Geom3DParser.h \
		mtxlib.h \
		RawivParser.h \
		reg3data.h \
		Reg3Parser.h

