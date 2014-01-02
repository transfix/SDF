// RawivPaser.cpp: implementation of the RawivPaser class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "RawivParser.h"

#include "reg3data.h"
#include "bio.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RawivParser::RawivParser()
{

}

RawivParser::~RawivParser()
{

}

bool RawivParser::parse(Reg3Data<float>* data, const char* fname)
{
	int nverts, ncells;
	float minext[3], maxext[3];
	
	if (!isRawivFile(fname)) {
		fprintf(stderr, "%s is not a rawiv file\n", fname);
		return false;
	}
	// determine file data type
	struct stat filestat;
	if (stat(fname, &filestat) < 0) {
		fprintf(stderr, "cannot find data file\n");
		return false;
	}
	int sz = filestat.st_size; 
	
	DiskIO *pio = new BufferedIO(fname);
	if(!pio->open()) {
		error("Data File Open Failed");
	}

	pio->get(minext, 3);
	pio->get(maxext, 3);
	pio->get(&nverts, 1);
	pio->get(&ncells, 1);
	pio->get(data->m_dim, 3);
	pio->get(data->m_orig, 3);
	pio->get(data->m_span, 3);

	nverts = (data->m_dim[0])*(data->m_dim[1])*(data->m_dim[2]); 
	if(data->p_data) { delete[] data->p_data; }
	data->p_data = new float[nverts];
	pio->get(data->p_data, nverts);
	
	data->init();
#ifdef _DEBUG
	printf("dim: %d %d %d\n", data->m_dim[0], data->m_dim[1], data->m_dim[2]);
	printf("orig: %f %f %f\n", data->m_orig[0], data->m_orig[1], data->m_orig[2]);
	printf("span: %f %f %f\n", data->m_span[0], data->m_span[1], data->m_span[2]);
#endif

	pio->close();
	delete pio;

	return true;
}

bool RawivParser::write(const Reg3Data<float>& data, const char* fname)
{
	DiskIO *pio = new BufferedIO(fname, DiskIO::WRITE);
	if(!pio->open()) {
		error("Cannot Open Data File to write");
		return false;
	}
	float maxext[3];
	maxext[0] = data.m_orig[0] + (data.m_dim[0]-1)*data.m_span[0];
	maxext[1] = data.m_orig[1] + (data.m_dim[1]-1)*data.m_span[1];
	maxext[2] = data.m_orig[2] + (data.m_dim[2]-1)*data.m_span[2];
	int nverts = data.getNVerts();
	int ncells = data.getNCells();
	pio->put(data.m_orig, 3);
	pio->put(maxext, 3);
	pio->put(&nverts, 1);
	pio->put(&ncells, 1);
	pio->put(data.m_dim, 3);
	pio->put(data.m_orig, 3);
	pio->put(data.m_span, 3);

	pio->put(data.p_data, nverts);

	pio->close(false);
	return true;
}
