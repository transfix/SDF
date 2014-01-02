

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


// SceneArrayNode.cpp: implementation of the SceneArrayNode class.
//
//////////////////////////////////////////////////////////////////////

#include "SceneArrayNode.h"
#include "Quaternion.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SceneArrayNode::SceneArrayNode(GLuint initializationCallList, GLuint finalizationCallList)
{
	m_Geometry = new Geometry;
	m_InitializationCallList = initializationCallList;
	m_FinalizationCallList = finalizationCallList;
}

SceneArrayNode::SceneArrayNode( Geometry* geometry, GLuint initializationCallList, GLuint finalizationCallList )
{
	this->m_Geometry = geometry;
	this->m_InitializationCallList = initializationCallList;
	this->m_FinalizationCallList = finalizationCallList;
}

SceneArrayNode::~SceneArrayNode()
{
	delete m_Geometry;
}

void SceneArrayNode::translate( float tx, float ty, float tz )
{
	m_Transformation.preMultiplication(CCVOpenGLMath::Matrix::translation(tx,ty,tz));
}

void SceneArrayNode::rotation( float angle, float x, float y, float z )
{
	m_Transformation.preMultiplication(CCVOpenGLMath::Quaternion::rotation(angle, x,y,z).buildMatrix());
}

void SceneArrayNode::scale( float s )
{
	m_Transformation.preMultiplication(CCVOpenGLMath::Matrix::scale(s,s,s));
}


void SceneArrayNode::setGeometry( Geometry* geometry )
{
	delete m_Geometry;
	m_Geometry = geometry;
}

Geometry* SceneArrayNode::getGeometry()
{
	return m_Geometry;
}

const Geometry* SceneArrayNode::getGeometry() const
{
	return m_Geometry;
}


void SceneArrayNode::setTransformation( const CCVOpenGLMath::Matrix& transformation)
{
	m_Transformation = transformation;
}

CCVOpenGLMath::Matrix& SceneArrayNode::getTransformation()
{
	return m_Transformation;
}

const CCVOpenGLMath::Matrix& SceneArrayNode::getTransformation() const
{
	return m_Transformation;
}
	
void SceneArrayNode::setInitializationCallList(GLuint initializationCallList)
{
	m_InitializationCallList = initializationCallList;
}

GLuint SceneArrayNode::getInitializationCallList() const
{
	return m_InitializationCallList;
}


void SceneArrayNode::setFinalizationCallList(GLuint finalizationCallList)
{
	m_FinalizationCallList = finalizationCallList;
}

GLuint SceneArrayNode::getFinalizationCallList() const
{
	return m_FinalizationCallList;
}

