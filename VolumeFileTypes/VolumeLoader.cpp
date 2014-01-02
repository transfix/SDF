/*****************************************************************************/
/*                             ______________________                        */
/*                            / _ _ _ _ _ _ _ _ _ _ _)                       */
/*            ____  ____  _  / /__  __  _____  __                            */
/*           (_  _)( ___)( \/ /(  \/  )(  _  )(  )                           */
/*             )(   )__)  )  (  )    (  )(_)(  )(__                          */
/*            (__) (____)/ /\_)(_/\/\_)(_____)(____)                         */
/*            _ _ _ _ __/ /                                                  */
/*           (___________/                     ___  ___                      */
/*                                      \  )| |   ) _ _|\   )                */
/*                                 ---   \/ | |  / |___| \_/                 */
/*                                                       _/                  */
/*                                                                           */
/*   Copyright (C) The University of Texas at Austin                         */
/*                                                                           */
/*     Author:     Vinay Siddavanahalli <skvinay@cs.utexas.edu>   2004-2005  */
/*                                                                           */
/*     Principal Investigator: Chandrajit Bajaj <bajaj@ices.utexas.edu>      */
/*                                                                           */
/*         Professor of Computer Sciences,                                   */
/*         Computational and Applied Mathematics Chair in Visualization,     */
/*         Director, Computational Visualization Center (CVC),               */
/*         Institute of Computational Engineering and Sciences (ICES)        */
/*         The University of Texas at Austin,                                */
/*         201 East 24th Street, ACES 2.324A,                                */
/*         1 University Station, C0200                                       */
/*         Austin, TX 78712-0027                                             */
/*         http://www.cs.utexas.edu/~bajaj                                   */
/*                                                                           */
/*         http://www.ices.utexas.edu/CVC                                    */
/*  This software comes with a license. Using this code implies that you     */
/*  read, understood and agreed to all the terms and conditions in that      */
/*  license.                                                                 */
/*                                                                           */
/*  We request that you agree to acknowledge the use of the software that    */
/*  results in any published work, including scientific papers, films and    */
/*  videotapes by citing the reference listed below                          */
/*                                                                           */
/*    C. Bajaj, P. Djeu, V. Siddavanahalli, A. Thane,                        */
/*    Interactive Visual Exploration of Large Flexible Multi-component       */
/*    Molecular Complexes,                                                   */
/*    Proc. of the Annual IEEE Visualization Conference, October 2004,       */
/*    Austin, Texas, IEEE Computer Society Press, pp. 243-250.               */
/*                                                                           */
/*****************************************************************************/// VolumeLoader.cpp: implementation of the VolumeLoader class.
//
//////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "VolumeLoader.h"

#include "RawIVFile.h"
#include "RawVFile.h"
#include "DXFile.h"
#include "MRCFile.h"
#include "VolumeFileType.h"
#include <string>
#include <cstring>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

VolumeLoader::VolumeLoader()
{
	// add each volume type to the maps
	addVolumeFileType(RawIVFile::getRepresentative());
	addVolumeFileType(RawVFile::getRepresentative());
	addVolumeFileType(DXFile::getRepresentative());
	addVolumeFileType(MRCFile::getRepresentative());
}

VolumeLoader::~VolumeLoader()
{

}

bool VolumeLoader::endsWith(string str, string substr)
{
	if( str.length() < substr.length() ) return false;
	if( str.substr(str.length() - substr.length(), substr.length()) == substr ) return true;
	return false;
}

bool VolumeLoader::saveFile(const string& fileName, SimpleVolumeData* simpleVolumeData)
{
	if( endsWith( fileName, "rawiv" ) )
		return saveFile( fileName, "RawIV files (*.rawiv)", simpleVolumeData);
	if( endsWith( fileName, "rawv" ) )
		return saveFile( fileName, "RawV files (*.rawv)", simpleVolumeData);
	if( endsWith( fileName, "dx" ) )
		return saveFile( fileName, "DX files (*.dx)", simpleVolumeData);
	if( endsWith( fileName, "mrc" ) )
		return saveFile( fileName, "MRC files (*.mrc)", simpleVolumeData);
	// try simple one
	return saveFile( fileName, "RawIV files (*.rawiv)", simpleVolumeData);
}

bool VolumeLoader::saveFile(const string& fileName, const string& selectedFilter, SimpleVolumeData* simpleVolumeData, unsigned int variable)
{
	if (m_FilterMap.count(selectedFilter) != 0) {
		string extension;
		string longName;

		if (fileName.rfind('.') != string::npos)
			extension = fileName.substr(fileName.rfind('.')+1);

		VolumeFileType* type = m_FilterMap[selectedFilter];

		// if no extension, add one
		if (extension.empty()) {
			longName = fileName + "." + type->extension();
		}
		else {
			longName = fileName;
		}

		return type->saveFile(simpleVolumeData, longName, variable);
	} 
	else {
		return false;
	}
}

SimpleVolumeData* VolumeLoader::loadFile(const string& fileName)
{
	string extension;

	if (fileName.rfind('.') != string::npos)
		extension = fileName.substr(fileName.rfind('.')+1);

	if (extension.empty() || m_ExtensionMap.count(extension) == 0) {
		// test every file type to find the correct one
		return tryAll(fileName);
	}
	else {
		// try to load the file with the correct file type
		VolumeFileType* type = m_ExtensionMap[extension];
		SimpleVolumeData* simpleVolumeData = type->loadFile(fileName);
		if (!simpleVolumeData) { // failed, try the other loaders
			return tryAll(fileName);
		}
		else {
			// success
			return simpleVolumeData;
		}
	}
}

string VolumeLoader::getLoadFilterString()
{
	string str("All Geometry Files ");
	str.append(getAllExtensions());
	// iterate through each loader and combine all the filters
	std::map<string, VolumeFileType*>::iterator it;
	for (it = m_FilterMap.begin(); it!=m_FilterMap.end(); ++it) {
		str.append(";;" + (it->first));
	}

	return str;
}

string VolumeLoader::getSaveFilterString()
{
	string str;
	bool first = true;
	// iterate through each loader and combine all the filters
	std::map<string, VolumeFileType*>::iterator it;
	for (it = m_FilterMap.begin(); it!=m_FilterMap.end(); ++it) {
		if (first) {
			first = false;
			str = it->first;
		}
		else 
			str.append(";;" + (it->first));
	}

	return str;
}

string VolumeLoader::getAllExtensions()
{
	string str("(");
	bool first = true;
	// iterate through each loader and combine all the filters
	std::map<string, VolumeFileType*>::iterator it;
	for (it = m_FilterMap.begin(); it!=m_FilterMap.end(); ++it) {
		if (first) {
			first = false;
			str.append("*." + (it->first));
		}
		else 
			str.append(" *." + (it->first));
	}
	str.append(")");

	return str;
}

SimpleVolumeData* VolumeLoader::tryAll(const string& fileName)
{
	// iterate through each loader and call checkType to determine
	// which loader can load the file
	std::map<string, VolumeFileType*>::iterator it;
	for (it = m_ExtensionMap.begin(); it!=m_ExtensionMap.end(); ++it) {
		if ((it->second)->checkType(fileName)) { // found it
			return (it->second)->loadFile(fileName);
		}
	}

	// didnt find the right loader
	return 0;
}

void VolumeLoader::addVolumeFileType(VolumeFileType* type)
{
	m_ExtensionMap[type->extension()] = type;
	m_FilterMap[type->filter()] = type;
}

bool VolumeLoader::isValidExtension( string extension )
{
	std::map<string, VolumeFileType*>::iterator it;
	for (it = m_ExtensionMap.begin(); it!=m_ExtensionMap.end(); ++it) {
		if( strcmp( (it->first).c_str(), extension.c_str()) == 0 )
			return true;
	}
	return false;
}

