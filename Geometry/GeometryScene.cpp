

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


// GeometryScene.cpp: implementation of the GeometryScene class.
//
//////////////////////////////////////////////////////////////////////

#include "GeometryScene.h"
#include "Geometry.h"

const unsigned int InitialArraySize = 16;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
GeometryScene::GeometryScene()
{
	initArrays();
	m_ReferenceCount = 1;
}

GeometryScene::~GeometryScene()
{
	deleteArrays();
}

void GeometryScene::initArrays() {

	initObjectArray();
	initIndexArray();
}

void GeometryScene::initObjectArray()
{
	m_NumberOfObjects = 0;
	m_SizeOfObjectsArray = InitialArraySize;

	m_SceneArrayNode = new SceneArrayNode*[m_SizeOfObjectsArray];
	m_ObjToIndex = new int[m_SizeOfObjectsArray];

	unsigned int c;
	for(c=0;c<m_SizeOfObjectsArray;c++) {
		m_SceneArrayNode[c] = 0;
	}
}

void GeometryScene::initIndexArray()
{
	m_NextIndexEntry = 0;
	m_SizeOfIndexToObjectsArray = InitialArraySize;
	m_IndexToObj = new int[m_SizeOfObjectsArray];

}

void GeometryScene::deleteArrays()
{
	delete [] m_SceneArrayNode;
	delete [] m_IndexToObj;
	delete [] m_ObjToIndex;
}

void GeometryScene::doubleObjectArray()
{
	unsigned int c;
	if (m_NumberOfObjects >= m_SizeOfObjectsArray) {
		SceneArrayNode** oldSceneArrayNode = m_SceneArrayNode;
		int* oldObjToIndex = m_ObjToIndex;
		m_SceneArrayNode = new SceneArrayNode*[m_SizeOfObjectsArray*2];
		m_ObjToIndex = new int[m_SizeOfObjectsArray*2];

		for(c=0;c<m_SizeOfObjectsArray;c++) {
			m_SceneArrayNode[c] = oldSceneArrayNode[c];
			m_ObjToIndex[c] = oldObjToIndex[c];
		}
		m_SizeOfObjectsArray *= 2;

		delete [] oldSceneArrayNode;
		delete [] oldObjToIndex;
	}
}

void GeometryScene::doubleIndexArray()
{
	unsigned int c;
	if (m_NextIndexEntry >= m_SizeOfIndexToObjectsArray) {
		int* oldIndexToObjectsArray = m_IndexToObj;
	
		m_IndexToObj = new int[m_SizeOfIndexToObjectsArray*2];

		for(c=0;c<m_NextIndexEntry;c++) {
			m_IndexToObj[c] = oldIndexToObjectsArray[c];
		}
		m_SizeOfIndexToObjectsArray *= 2;

		delete [] oldIndexToObjectsArray;
	}
}

int GeometryScene::add( SceneArrayNode* sceneArrayNode )
{
	doubleObjectArray();
	m_SceneArrayNode[m_NumberOfObjects] = sceneArrayNode;
	int indexPosition;
	if (m_HoleList.isEmpty()) {
		doubleIndexArray();
		indexPosition = m_NextIndexEntry;
		m_NextIndexEntry++;
	}
	else {
		indexPosition = m_HoleList.deQueue();
	}
	m_ObjToIndex[m_NumberOfObjects] = indexPosition;
	m_IndexToObj[indexPosition] = m_NumberOfObjects;
	m_NumberOfObjects++;
	return m_NumberOfObjects-1;
}

int GeometryScene::add( Geometry* geometry )
{
	SceneArrayNode* sceneArrayNode = new SceneArrayNode( geometry );
	return add(sceneArrayNode);
}

bool GeometryScene::set( SceneArrayNode* sceneArrayNode, unsigned int index )
{
	index = m_IndexToObj[index];
	if( index < m_NumberOfObjects && index >=0) 
	{
		m_SceneArrayNode[index] = sceneArrayNode;
		return true;
	}
	return false;	
}

SceneArrayNode* GeometryScene::remove( unsigned int index )
{
	int object = m_IndexToObj[index];
	if( object < (int)m_NumberOfObjects && object >=0 ) 
	{
		m_HoleList.enQueue(index);
		m_NumberOfObjects--;
		SceneArrayNode* temp = m_SceneArrayNode[object];
		m_SceneArrayNode[object] = m_SceneArrayNode[m_NumberOfObjects];
		m_ObjToIndex[object] = m_ObjToIndex[m_NumberOfObjects];
		m_IndexToObj[m_ObjToIndex[object]] = object;
		return temp;
	}
	else
	{
		return 0;
	}
}

SceneArrayNode* GeometryScene::get( unsigned int index )
{
	int object = m_IndexToObj[index];
	if( object < (int)m_NumberOfObjects && object >=0 ) 
	{
		return m_SceneArrayNode[object];
	}
	else
	{
		return 0;
	}
}

SceneArrayNode* GeometryScene::getIth(unsigned int i) const
{
	if( i < m_NumberOfObjects && i >=0 ) 
	{
		return m_SceneArrayNode[i];
	}
	else
	{
		return 0;
	}
}

unsigned int GeometryScene::getNumberOfSceneArrayNodes()
{
	return m_NumberOfObjects;
}

void GeometryScene::clear()
{
	unsigned int c;
	for(c=0;c<m_NumberOfObjects;c++)
	{
		delete m_SceneArrayNode[c];
		m_SceneArrayNode[c] = 0;
	}
	m_NumberOfObjects=0;
	m_NextIndexEntry = 0;
	m_HoleList.clearQueue();
}

void GeometryScene::addMe()
{
	m_ReferenceCount++;
}

void GeometryScene::deleteMe()
{
	m_ReferenceCount--;
	if( m_ReferenceCount == 0 ) 
	{
		clear();
		delete this;
	}
}

