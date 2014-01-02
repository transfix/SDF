// RawcFile.h: interface for the RawcFile class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RAWCFILE_H__F47E5524_8A5C_43D2_AFD5_8E6E5ACCB23C__INCLUDED_)
#define AFX_RAWCFILE_H__F47E5524_8A5C_43D2_AFD5_8E6E5ACCB23C__INCLUDED_

#include "GeometryFileType.h"

class RawcFile : public GeometryFileType  
{
public:
	virtual ~RawcFile();

	virtual Geometry* loadFile(const string& fileName);
	virtual bool checkType(const string& fileName);
	virtual bool saveFile(const Geometry* geometry, const string& fileName);

	virtual string extension() { return "rawc"; };
	virtual string filter() { return "Rawc files (*.rawc)"; };

	static RawcFile ms_RawcFileRepresentative;
	static GeometryFileType* getRepresentative();

protected:
	RawcFile();

};

#endif // !defined(AFX_RAWCFILE_H__F47E5524_8A5C_43D2_AFD5_8E6E5ACCB23C__INCLUDED_)
