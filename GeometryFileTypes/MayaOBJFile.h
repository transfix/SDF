// MayaOBJFile.h: interface for the MayaOBJFile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MAYAOBJFILE_H__1A1D22F0_B170_4AD3_8BDA_35CF0F0C54B4__INCLUDED_)
#define AFX_MAYAOBJFILE_H__1A1D22F0_B170_4AD3_8BDA_35CF0F0C54B4__INCLUDED_

#include "GeometryFileType.h"

class MayaOBJFile : public GeometryFileType  
{
public:
	virtual ~MayaOBJFile();

	virtual Geometry* loadFile(const string& fileName);
	virtual bool checkType(const string& fileName);
	virtual bool saveFile(const Geometry* geometry, const string& fileName);

	virtual string extension() { return "obj"; };
	virtual string filter() { return "Maya OBJ files (*.obj)"; };

	static MayaOBJFile ms_MayaOBJFileRepresentative;
	static GeometryFileType* getRepresentative();

protected:
	MayaOBJFile();
};

#endif // !defined(AFX_MAYAOBJFILE_H__1A1D22F0_B170_4AD3_8BDA_35CF0F0C54B4__INCLUDED_)
