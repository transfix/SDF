/*
  Copyright 2002-2003 The University of Texas at Austin
  
	Authors: Anthony Thane <thanea@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with iotree; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

// GeometryLoader.cpp: implementation of the GeometryLoader class.
//
//////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "GeometryLoader.h"
#include "RawFile.h"
#include "RawnFile.h"
#include "RawcFile.h"
#include "RawncFile.h"
#include "MayaOBJFile.h"
//#include "C2CFile.h"
#include "GeometryFileType.h"
#include <string>
#include <cstring>
//#include <qstring.h>
//#include <qfileinfo.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GeometryLoader::GeometryLoader()
{
	// add each geometry type to the maps
//	addGeometryFileType(C2CFile::getRepresentative());
	addGeometryFileType(MayaOBJFile::getRepresentative());
	addGeometryFileType(RawncFile::getRepresentative());
	addGeometryFileType(RawcFile::getRepresentative());
	addGeometryFileType(RawnFile::getRepresentative());
	addGeometryFileType(RawFile::getRepresentative());

}

GeometryLoader::~GeometryLoader()
{

}

bool GeometryLoader::endsWith(string str, string substr)
{
	if( str.length() < substr.length() ) return false;
	if( str.substr(str.length() - substr.length(), substr.length()) == substr ) return true;
	return false;
}

bool GeometryLoader::saveFile(const string& fileName, Geometry* geometry)
{
	if( endsWith( fileName, "rawnc" ) )
		return saveFile( fileName, "Rawnc files (*.rawnc)", geometry);
	if( endsWith( fileName, "rawn" ) )
		return saveFile( fileName, "Rawn files (*.rawn)", geometry);
	if( endsWith( fileName, "rawc" ) )
		return saveFile( fileName, "Rawc files (*.rawc)", geometry);
	if( endsWith( fileName, "raw" ) )
		return saveFile( fileName, "Raw files (*.raw)", geometry);
	if( endsWith( fileName, "c2c" ) )
		return saveFile( fileName, "c2c files (*.c2c)", geometry);
	if( endsWith( fileName, "obj" ) )
		return saveFile( fileName, "Maya OBJ files (*.obj)", geometry);

	return saveFile( fileName, "Rawnc files (*.rawnc)", geometry);
}

bool GeometryLoader::saveFile(const string& fileName, const string& selectedFilter, Geometry* geometry)
{
	//if (m_FilterMap.contains(selectedFilter)) {
	if (m_FilterMap.count(selectedFilter) != 0) {
		//QFileInfo fileInfo(fileName);
		//QString extension = fileInfo.extension(false);
		string extension;
		string longName;

		if (fileName.rfind('.') != string::npos)
			extension = fileName.substr(fileName.rfind('.')+1);

		GeometryFileType* type = m_FilterMap[selectedFilter];

		// if no extension, add one
		if (extension.empty()) {
			longName = fileName + "." + type->extension();
		}
		else {
			longName = fileName;
		}

		return type->saveFile(geometry, longName);
	} 
	else {
		return false;
	}
}

Geometry* GeometryLoader::loadFile(const string& fileName)
{
	//QFileInfo fileInfo(fileName);
	//QString extension = fileInfo.extension(false);
	string extension;

	if (fileName.rfind('.') != string::npos)
		extension = fileName.substr(fileName.rfind('.')+1);

	//printf("extension = %s\n", extension.c_str());
	
	//if (extension.isEmpty() || !m_ExtensionMap.contains(extension)) {
	if (extension.empty() || m_ExtensionMap.count(extension) == 0) {
		// test every file type to find the correct one
		return tryAll(fileName);
	}
	else {
		// try to load the file with the correct file type
		GeometryFileType* type = m_ExtensionMap[extension];
		Geometry* geometry = type->loadFile(fileName);
		if (!geometry) { // failed, try the other loaders
			return tryAll(fileName);
		}
		else {
			// success
			return geometry;
		}
	}
}

string GeometryLoader::getLoadFilterString()
{
	//QString string("All Geometry Files ");
	//string.append(getAllExtensions());
	string str("All Geometry Files ");
	str.append(getAllExtensions());
	// iterate through each loader and combine all the filters
	//QMap<QString, GeometryFileType*>::Iterator it;
	std::map<string, GeometryFileType*>::iterator it;
	for (it = m_FilterMap.begin(); it!=m_FilterMap.end(); ++it) {
		//string.append(";;" + it.key());
		str.append(";;" + (it->first));
	}

	return str;
}

string GeometryLoader::getSaveFilterString()
{
	//QString string;
	string str;
	bool first = true;
	// iterate through each loader and combine all the filters
	//QMap<QString, GeometryFileType*>::Iterator it;
	std::map<string, GeometryFileType*>::iterator it;
	for (it = m_FilterMap.begin(); it!=m_FilterMap.end(); ++it) {
		if (first) {
			first = false;
			//string = it.key();
			str = it->first;
		}
		else 
			//string.append(";;" + it.key());
			str.append(";;" + (it->first));
	}

	return str;
}

string GeometryLoader::getAllExtensions()
{
	//QString string("(");
	string str("(");
	bool first = true;
	// iterate through each loader and combine all the filters
	//QMap<QString, GeometryFileType*>::Iterator it;
	std::map<string, GeometryFileType*>::iterator it;
	for (it = m_FilterMap.begin(); it!=m_FilterMap.end(); ++it) {
		if (first) {
			first = false;
			//string.append("*." + it.key());
			str.append("*." + (it->first));
		}
		else 
			//string.append(" *." + it.key());
			str.append(" *." + (it->first));
	}
	//string.append(")");
	str.append(")");

	return str;
}

Geometry* GeometryLoader::tryAll(const string& fileName)
{
	// iterate through each loader and call checkType to determine
	// which loader can load the file
	//QMap<QString, GeometryFileType*>::Iterator it;
	std::map<string, GeometryFileType*>::iterator it;
	for (it = m_ExtensionMap.begin(); it!=m_ExtensionMap.end(); ++it) {
		//if (it.data()->checkType(fileName)) { // found it
		//	return it.data()->loadFile(fileName);
		//}
		if ((it->second)->checkType(fileName)) { // found it
			return (it->second)->loadFile(fileName);
		}
	}

	// didnt find the right loader
	return 0;
}

void GeometryLoader::addGeometryFileType(GeometryFileType* type)
{
	m_ExtensionMap[type->extension()] = type;
	m_FilterMap[type->filter()] = type;
}

bool GeometryLoader::isValidExtension( string extension )
{
	std::map<string, GeometryFileType*>::iterator it;
	for (it = m_ExtensionMap.begin(); it!=m_ExtensionMap.end(); ++it) {
		if( strcmp( (it->first).c_str(), extension.c_str()) == 0 )
			return true;
	}
	return false;
}
