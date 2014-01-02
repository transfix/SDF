/*
  Copyright (c): Xiaoyu Zhang (xiaoyu@csusm.edu)

  This file is part of sdf (signed distance function).

  sdf is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  sdf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with sdf; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
// RawivPaser.h: interface for the RawivPaser class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RAWIVPASER_H__938E0E03_F02E_4298_B690_AD4DE7EF3AB4__INCLUDED_)
#define AFX_RAWIVPASER_H__938E0E03_F02E_4298_B690_AD4DE7EF3AB4__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string.h>

#include "Reg3Parser.h"

class RawivParser : public Reg3Parser  
{
public:
	RawivParser();
	virtual ~RawivParser();

	virtual bool parse(Reg3Data<float>* data, const char* fname);

	virtual bool write(const Reg3Data<float>& data, const char* fname);

private:
	bool isRawivFile(const char* fname) {
		int len = (int)strlen(fname);
		return (len > 6 && strcmp(fname+len-6, ".rawiv") == 0);
	}
};

#endif // !defined(AFX_RAWIVPASER_H__938E0E03_F02E_4298_B690_AD4DE7EF3AB4__INCLUDED_)
