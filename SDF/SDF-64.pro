TEMPLATE = app
CONFIG  += opengl warn_off
#TARGET  += sdf 

macx-g++ {
  QMAKE_CFLAGS += -m64
  QMAKE_CXXFLAGS += -m64
}


INCLUDEPATH += ../ByteOrder ../GeometryFileTypes ../UsefulMath ../Geometry ../SignDistanceFunction ../VolumeFileTypes

LIBS += ../ByteOrder/libByteOrder.a ../GeometryFileTypes/libGeometryFileTypes.a ../UsefulMath/libUsefulMath.a ../Geometry/libGeometry.a ../SignDistanceFunction/libSignDistanceFunction.a ../VolumeFileTypes/libVolumeFileTypes.a

# Input

SOURCES = main.cpp	
