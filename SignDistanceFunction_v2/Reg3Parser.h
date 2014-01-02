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
// Reg3Parser.h: interface for the Reg3Parser class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_REG3PARSER_H__AE88C884_6FE3_433E_B6D2_F84CD6E39724__INCLUDED_)
#define AFX_REG3PARSER_H__AE88C884_6FE3_433E_B6D2_F84CD6E39724__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "reg3data.h"

class Reg3Parser  
{
public:
	Reg3Parser();
	virtual ~Reg3Parser();

	virtual bool parse(Reg3Data<float>* data, const char* fname) {
		return true;
	}

	virtual bool write(const Reg3Data<float>& data, const char* fname) {
		return true;
	}
};

#endif // !defined(AFX_REG3PARSER_H__AE88C884_6FE3_433E_B6D2_F84CD6E39724__INCLUDED_)
