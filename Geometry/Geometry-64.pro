# qmake project generated by QMsDev
#
# General settings

TEMPLATE = lib
CONFIG  += opengl warn_off staticlib create_prl
#TARGET  += Geometry

macx-g++ {
  QMAKE_CFLAGS += -m64
  QMAKE_CXXFLAGS += -m64
}


# Input
INCLUDEPATH += ../UsefulMath

SOURCES =  \
		Geometry.cpp \
		GeometryScene.cpp \
		GeometrySceneArray.cpp \
		IntQueue.cpp \
		MeshDerivatives.cpp \
		SceneArrayNode.cpp \
		Texture2D.cpp

HEADERS =  \
		Geometry.h \
		GeometryScene.h \
		GeometrySceneArray.h \
		IntQueue.h \
		MeshDerivatives.h \
		SceneArrayNode.h \
		Texture2D.h

