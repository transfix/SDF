

/***********************************************************************************/
/*																				   */
/*	  Copyright 2003 University of Texas at Austin                                 */
/*	  Supervisor: Dr C Bajaj bajaj@cs.utexas.edu,                                  */
/*    Authors:    Anthony Thane thanea@ices.utexas.edu                             */
/*                S K Vinay  skvinay@cs.utexas.edu                                 */
/*																				   */
/*    This program is free software; you can redistribute it and/or modify         */
/*    it under the terms of the GNU General Public License as published by         */
/*    the Free Software Foundation; either version 2 of the License, or            */
/*    (at your option) any later version.                                          */
/*																				   */
/*    This program is distributed in the hope that it will be useful,              */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of               */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                */
/*    GNU General Public License for more details.                                 */
/*																				   */
/*    You should have received a copy of the GNU General Public License			   */
/*    along with this program; if not, write to the Free Software                  */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    */
/*                                                                                 */
/***********************************************************************************/


// GeometrySceneArray.h: interface for the GeometrySceneArray class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GEOMETRYSCENEARRAY_H__A5F41949_2CB3_476B_9C01_F48B9C11E564__INCLUDED_)
#define AFX_GEOMETRYSCENEARRAY_H__A5F41949_2CB3_476B_9C01_F48B9C11E564__INCLUDED_

class GeometryScene;

class GeometrySceneArray
{
public:
	GeometrySceneArray();
	virtual ~GeometrySceneArray();

	void initArray();

	void doubleArray();
	int add( GeometryScene* geometryScene );
	bool set( GeometryScene* geometryScene, unsigned int index );
	GeometryScene* remove( unsigned int index );
	GeometryScene* get( unsigned int index );
	unsigned int getNumberOfObjects();

	void clear();

protected:
	GeometryScene** m_GeometriesScene;
	unsigned int m_NumberOfObjects;
	unsigned int m_SizeOfObjectsArray;
};


#endif // !defined(AFX_GEOMETRYSCENEARRAY_H__A5F41949_2CB3_476B_9C01_F48B9C11E564__INCLUDED_)
