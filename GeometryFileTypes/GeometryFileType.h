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

// GeometryFileType.h: interface for the GeometryFileType class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GEOMETRYFILETYPE_H__3D1E7183_20B6_4CBE_A5E4_F0C92642AB4E__INCLUDED_)
#define AFX_GEOMETRYFILETYPE_H__3D1E7183_20B6_4CBE_A5E4_F0C92642AB4E__INCLUDED_

//#include <qstring.h>
#include <string>
using std::string;

class Geometry;

class GeometryFileType  
{
public:
	GeometryFileType();
	virtual ~GeometryFileType();

	virtual Geometry* loadFile(const string& fileName) = 0;
	virtual bool checkType(const string& fileName) = 0;
	virtual bool saveFile(const Geometry* geometry, const string& fileName) = 0;

	virtual string extension() = 0;
	virtual string filter() = 0;

};

#endif // !defined(AFX_GEOMETRYFILETYPE_H__3D1E7183_20B6_4CBE_A5E4_F0C92642AB4E__INCLUDED_)
