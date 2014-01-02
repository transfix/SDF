TEMPLATE = app
CONFIG  += warn_off console
CONFIG -= qt
#TARGET  += sdf 

INCLUDEPATH += ../ByteOrder ../GeometryFileTypes ../UsefulMath ../Geometry ../SignDistanceFunction ../VolumeFileTypes

LIBS += ../ByteOrder/libByteOrder.a ../GeometryFileTypes/libGeometryFileTypes.a ../UsefulMath/libUsefulMath.a ../Geometry/libGeometry.a ../SignDistanceFunction/libSignDistanceFunction.a ../VolumeFileTypes/libVolumeFileTypes.a

# Input

SOURCES = main.cpp	
