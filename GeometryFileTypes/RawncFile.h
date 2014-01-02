// RawncFile.h: interface for the RawncFile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RAWNCFILE_H__2755A9D8_A22A_496D_BF36_C7CF6FE93F6D__INCLUDED_)
#define AFX_RAWNCFILE_H__2755A9D8_A22A_496D_BF36_C7CF6FE93F6D__INCLUDED_

#include "GeometryFileType.h"

class RawncFile : public GeometryFileType  
{
public:
	virtual ~RawncFile();

	virtual Geometry* loadFile(const string& fileName);
	virtual bool checkType(const string& fileName);
	virtual bool saveFile(const Geometry* geometry, const string& fileName);

	virtual string extension() { return "rawnc"; };
	virtual string filter() { return "Rawnc files (*.rawnc)"; };

	static RawncFile ms_RawncFileRepresentative;
	static GeometryFileType* getRepresentative();

	bool saveMolCenters(const Geometry* geometry);

protected:
	RawncFile();
};

#endif // !defined(AFX_RAWNCFILE_H__2755A9D8_A22A_496D_BF36_C7CF6FE93F6D__INCLUDED_)
