/*****************************************************************************/
/*                             ______________________                        */
/*                            / _ _ _ _ _ _ _ _ _ _ _)                       */
/*            ____  ____  _  / /__  __  _____  __                            */
/*           (_  _)( ___)( \/ /(  \/  )(  _  )(  )                           */
/*             )(   )__)  )  (  )    (  )(_)(  )(__                          */
/*            (__) (____)/ /\_)(_/\/\_)(_____)(____)                         */
/*            _ _ _ _ __/ /                                                  */
/*           (___________/                     ___  ___                      */
/*                                      \  )| |   ) _ _|\   )                */
/*                                 ---   \/ | |  / |___| \_/                 */
/*                                                       _/                  */
/*                                                                           */
/*   Copyright (C) The University of Texas at Austin                         */
/*                                                                           */
/*     Author:     Vinay Siddavanahalli <skvinay@cs.utexas.edu>   2004-2005  */
/*                                                                           */
/*     Principal Investigator: Chandrajit Bajaj <bajaj@ices.utexas.edu>      */
/*                                                                           */
/*         Professor of Computer Sciences,                                   */
/*         Computational and Applied Mathematics Chair in Visualization,     */
/*         Director, Computational Visualization Center (CVC),               */
/*         Institute of Computational Engineering and Sciences (ICES)        */
/*         The University of Texas at Austin,                                */
/*         201 East 24th Street, ACES 2.324A,                                */
/*         1 University Station, C0200                                       */
/*         Austin, TX 78712-0027                                             */
/*         http://www.cs.utexas.edu/~bajaj                                   */
/*                                                                           */
/*         http://www.ices.utexas.edu/CVC                                    */
/*  This software comes with a license. Using this code implies that you     */
/*  read, understood and agreed to all the terms and conditions in that      */
/*  license.                                                                 */
/*                                                                           */
/*  We request that you agree to acknowledge the use of the software that    */
/*  results in any published work, including scientific papers, films and    */
/*  videotapes by citing the reference listed below                          */
/*                                                                           */
/*    C. Bajaj, P. Djeu, V. Siddavanahalli, A. Thane,                        */
/*    Interactive Visual Exploration of Large Flexible Multi-component       */
/*    Molecular Complexes,                                                   */
/*    Proc. of the Annual IEEE Visualization Conference, October 2004,       */
/*    Austin, Texas, IEEE Computer Society Press, pp. 243-250.               */
/*                                                                           */
/*****************************************************************************/// RawVFile.cpp: implementation of the RawVFile class.
//
//////////////////////////////////////////////////////////////////////

#include "RawVFile.h"
#include "SimpleVolumeData.h"
#include "ByteSwapping.h"
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(_LARGEFILE_SOURCE)
# define FOPEN fopen64
#else
# define FOPEN fopen
#endif

RawVFile RawVFile::ms_RawVFileRepresentative;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RawVFile::RawVFile()
{

}

RawVFile::~RawVFile()
{

}

SimpleVolumeData* RawVFile::loadFile(const string& fileName)
{
	SimpleVolumeData* simpleVolumeData = new SimpleVolumeData(128, 128, 128);

	FILE *fp=FOPEN(fileName.c_str(), "rb");
	unsigned int i, magic, dims[3], numTimeSteps, numVariables;
	float minExt[4], maxExt[4];

	// check to make sure the file exists
	if (!fp) {
		printf("Error: could not open file\n");
		delete simpleVolumeData; simpleVolumeData = 0;
		return false;
	}

	// read the header
	//
	// read the magic number
	fread(&magic, 1, sizeof(unsigned int), fp);
	if (isLittleEndian()) swapByteOrder(&magic, 1);
	if (magic != 0xBAADBEEF) {
		// doh! this isn't a RawV file
		printf("Error: file format is not recognized. (expecting RawV)\n");
		fclose(fp);
		return false;
	}
	// the dimensions, # of timesteps, # of variables, and min/max extents
	fread(dims, 3, sizeof(unsigned int), fp);
	fread(&numTimeSteps, 1, sizeof(unsigned int), fp);
	fread(&numVariables, 1, sizeof(unsigned int), fp);
	fread(minExt, 4, sizeof(float), fp);
	fread(maxExt, 4, sizeof(float), fp);
	if (isLittleEndian()) {
		swapByteOrder(dims, 3);
		swapByteOrder(&numTimeSteps, 1);
		swapByteOrder(&numVariables, 1);
		swapByteOrder(minExt, 4);
		swapByteOrder(maxExt, 4);
	}

	// call a bunch of set...() functions
	//
	// setNumberOfVariables() has the convenient side effect of cleaning up
	// any data that this instance might already be holding. (effectively wiping
	// the slate clean)
	simpleVolumeData->setNumberOfVariables(numVariables);
	// the rest ...
	simpleVolumeData->setDimensions(dims);
	simpleVolumeData->setMinExtent(minExt);
	simpleVolumeData->setMaxExtent(maxExt);

	// finish reading the header (variable names and types)
	for (i=0; i < numVariables; i++) {
		char name[64];
		unsigned char type;

		// read the type and the name of a variable
		fread(&type, 1, sizeof(unsigned char), fp);
		fread(name, 64, sizeof(char), fp);
		name[63] = '\0'; // just in case
		// call set...() functions
		switch(type) {
			case 1:	simpleVolumeData->setType(i, SimpleVolumeData::UCHAR); break;
			case 2: simpleVolumeData->setType(i, SimpleVolumeData::USHORT); break;
			case 3: simpleVolumeData->setType(i, SimpleVolumeData::ULONG); break;
			case 4: simpleVolumeData->setType(i, SimpleVolumeData::FLOAT); break;
			case 5: simpleVolumeData->setType(i, SimpleVolumeData::DOUBLE); break;
			default: break;
		}
		simpleVolumeData->setName(i, name);
	}
	
	// finally, read the data
	for (i=0; i < numVariables; i++) {
		// now that all the meta-data is filled in, allocate space for the variables
		unsigned char* data = new unsigned char [simpleVolumeData->getDataSize(i)];
		// read the data for this variables first (and only?) timestep
		fread(data, simpleVolumeData->getDataSize(i), 1, fp);
		// swap the byte order if needed
		if (isLittleEndian()) {
			switch(simpleVolumeData->getTypeSize(i)) {
				case 1:
					break;
				case 2:
					swapByteOrder((unsigned short *)data, dims[0]*dims[1]*dims[2]);
					break;
				case 4:
					swapByteOrder((float *)data, dims[0]*dims[1]*dims[2]);
					break;
				case 8:
					swapByteOrder((double *)data, dims[0]*dims[1]*dims[2]);
					break;
				default:
					break;
			}
		}
		simpleVolumeData->setData(i, data);
		// seek past any remaining timesteps
		if (numTimeSteps > 1)
			fseek(fp, simpleVolumeData->getDataSize(i)*(numTimeSteps-1), SEEK_CUR);
	}
	
	// close the file
	fclose(fp);
	
	return simpleVolumeData;
}

bool RawVFile::checkType(const string& fileName)
{
	return false;
}

bool RawVFile::saveFile(SimpleVolumeData* simpleVolumeData, const string& fileName, unsigned int variable)
{
	// no data? bail.
	if (!simpleVolumeData->getData(variable))
		return false;

	// unsupported datatype for the requested variable? bail.
	if (simpleVolumeData->getType(variable) == SimpleVolumeData::ULONG || simpleVolumeData->getType(variable) == SimpleVolumeData::DOUBLE)
		return false;

	FILE *fp = FOPEN(fileName.c_str(), "wb");
	char *header=NULL;
	unsigned int hsize=0;
	
	// failed to open the file? bail.
	if (!fp)
		return false;

	// make the data big endian
	simpleVolumeData->makeVariablesBigEndian();

	// create the header
	simpleVolumeData->createRawVHeader(&header, &hsize);

	// write the header
	fwrite(header, hsize, 1, fp);

	// write the data
	for (unsigned int i=0; i < simpleVolumeData->getNumberOfVariables(); i++)
		fwrite(simpleVolumeData->getData(i), simpleVolumeData->getDataSize(i), 1, fp);

	// make the data native endian (poorly named function)
	simpleVolumeData->makeVariablesBigEndian();
	
	// clean up
	delete [] header;
	fclose(fp);
	
	return true;
}

VolumeFileType* RawVFile::getRepresentative()
{
	return &ms_RawVFileRepresentative;
}

