# qmake project generated by QMsDev
#
# General settings

TEMPLATE = lib
CONFIG  += warn_off staticlib create_prl
#TARGET  += VolumeFileTypes

# Input

INCLUDEPATH += ../ByteOrder
INCLUDEPATH += ../

SOURCES =  \
		DXFile.cpp \
		MRCFile.cpp \
		SimpleVolumeData.cpp \
		RawIVFile.cpp \
		RawVFile.cpp \
		VolumeFileType.cpp \
		VolumeLoader.cpp

HEADERS =  \
		DXFile.h \
		MRCFile.h \
		SimpleVolumeData.h \
		RawIVFile.h \
		RawVFile.h \
		VolumeFileType.h \
		VolumeLoader.h 

