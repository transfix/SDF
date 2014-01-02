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

// RawFile.h: interface for the RawFile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RAWFILE_H__C647E6C6_7669_4CFE_AFBC_BB072A6A8EC7__INCLUDED_)
#define AFX_RAWFILE_H__C647E6C6_7669_4CFE_AFBC_BB072A6A8EC7__INCLUDED_

#include "GeometryFileType.h"

//using std::string;

class RawFile : public GeometryFileType  
{
public:
	virtual ~RawFile();

	virtual Geometry* loadFile(const string& fileName);
	virtual bool checkType(const string& fileName);
	virtual bool saveFile(const Geometry* geometry, const string& fileName);

	virtual string extension() { return "raw"; };
	virtual string filter() { return "Raw files (*.raw)"; };

	static RawFile ms_RawFileRepresentative;
	static GeometryFileType* getRepresentative();

protected:
	RawFile();


};

#endif // !defined(AFX_RAWFILE_H__C647E6C6_7669_4CFE_AFBC_BB072A6A8EC7__INCLUDED_)
