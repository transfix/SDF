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
/*****************************************************************************/// RawIVFile.cpp: implementation of the RawIVFile class.
//
//////////////////////////////////////////////////////////////////////

#include "RawIVFile.h"
#include "SimpleVolumeData.h"
#include "ByteSwapping.h"
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdlib>

#if defined(_LARGEFILE_SOURCE)
# define FOPEN fopen64
#else
# define FOPEN fopen
#endif

RawIVFile RawIVFile::ms_RawIVFileRepresentative;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RawIVFile::RawIVFile()
{

}

RawIVFile::~RawIVFile()
{

}

SimpleVolumeData* RawIVFile::loadFile(const string& fileName)
{
	SimpleVolumeData* simpleVolumeData = new SimpleVolumeData(128, 128, 128);

	FILE *fp=FOPEN(fileName.c_str(), "rb");
	unsigned int dims[3];
	unsigned int temp_numverts, temp_numcells; //for large files, these are meaningless, and anyway not used.
	float minExt[3], maxExt[3], orig[3], span[3];

	// check to make sure the file exists
	if (!fp) {
		printf("Error: could not open file\n");
		delete simpleVolumeData; simpleVolumeData = 0;
		return 0;
	}

	// read the header
	// 
	fread(minExt, 3, sizeof(float), fp);
	fread(maxExt, 3, sizeof(float), fp);
	fread(&temp_numverts, 1, sizeof(unsigned int), fp);
	fread(&temp_numcells, 1, sizeof(unsigned int), fp);
	fread(dims, 3, sizeof(unsigned int), fp);
	fread(orig, 3, sizeof(float), fp);
	fread(span, 3, sizeof(float), fp);
	if (isLittleEndian()) {
		swapByteOrder(minExt, 3);
		swapByteOrder(maxExt, 3);
		swapByteOrder(&temp_numverts, 1);
		swapByteOrder(&temp_numcells, 1);
		swapByteOrder(dims, 3);
		swapByteOrder(orig, 3);
		swapByteOrder(span, 3);
	}

	/*	//// only for NMJ project! //////
		{
			minExt[0] *= 150;
			minExt[1] *= 150;
			minExt[2] *= 150;
		}
		{
			maxExt[0] *= 150;
			maxExt[1] *= 150;
			maxExt[2] *= 150;

			orig[0] *= 150;
			orig[1] *= 150;
			orig[2] *= 150;

			span[0] = span[1] = span[2] = 150;
		}
	*/	/////////////////////////////////


	// find out how large the data is
	Q_ULLONG dataStart = ftell(fp), dataSize;
	fseek(fp, 0, SEEK_END);
	dataSize = ftell(fp) - dataStart;
	fseek(fp, dataStart, SEEK_SET);

	// a small sanity check to make sure this file is valid
	if ((dataSize % ((Q_ULLONG)dims[0]*dims[1]*dims[2])) != 0) {
		printf("Error: rawiv file header dimensions don't match file size");
		fclose(fp);
		return false;
	}
	
	// call some set...() functions
	//
	simpleVolumeData->setNumberOfVariables(1);
	simpleVolumeData->setDimensions(dims);
	simpleVolumeData->setMinExtent(minExt);
	simpleVolumeData->setMaxExtent(maxExt);

	// figure out the data type
	switch (dataSize / ((Q_ULLONG)dims[0]*dims[1]*dims[2])) 
	{
	case 1: simpleVolumeData->setType(0, SimpleVolumeData::UCHAR); break;
	case 2: simpleVolumeData->setType(0, SimpleVolumeData::USHORT); break;
	case 4: simpleVolumeData->setType(0, SimpleVolumeData::FLOAT); break;
	default: simpleVolumeData->setType(0, SimpleVolumeData::NO_TYPE); break;
	}

	// allocate space for the data
	unsigned char * data = (unsigned char*)malloc((Q_ULLONG)sizeof(unsigned char)*dataSize);

	// read the data
	fread(data, dataSize, 1, fp);
	// swap the byte order if needed
	if (isLittleEndian()) {
		switch(simpleVolumeData->getTypeSize(0)) {
			case 1:
				break;
			case 2:
				swapByteOrder((unsigned short *)data, (Q_ULLONG)dims[0]*dims[1]*dims[2]);
				break;
			case 4:
				swapByteOrder((float *)data, (Q_ULLONG)dims[0]*dims[1]*dims[2]);
				break;
			case 8:
				swapByteOrder((double *)data, (Q_ULLONG)dims[0]*dims[1]*dims[2]);
				break;
			default:
				break;
		}
	}

	simpleVolumeData->setData(0, data);
	// close the file

	fclose(fp);
	return simpleVolumeData;
}

bool RawIVFile::checkType(const string& fileName)
{
	return false;
}

bool RawIVFile::saveFile(SimpleVolumeData* simpleVolumeData, const string& fileName, unsigned int variable)
{
	// no data? bail.
	if (!simpleVolumeData->getData(variable))
		return false;

	// unsupported datatype for the requested variable? bail.
	if (simpleVolumeData->getType(variable) == SimpleVolumeData::ULONG || simpleVolumeData->getType(variable) == SimpleVolumeData::DOUBLE)
		return false;

	FILE *fp = FOPEN(fileName.c_str(), "wb");
	char header[68];
	
	// failed to open the file? bail.
	if (!fp)
		return false;

	// make the data big endian
	simpleVolumeData->makeVariablesBigEndian();

	// create the header
	simpleVolumeData->createRawIVHeader(header);

	// write the header
	fwrite(header, sizeof(header), 1, fp);


	// write the data
	fwrite(simpleVolumeData->getData(variable), simpleVolumeData->getDataSize(variable), 1, fp);

	// make the data native endian (poorly named function)
	simpleVolumeData->makeVariablesBigEndian();
	
	// clean up
	fclose(fp);
	
	return true;
}

VolumeFileType* RawIVFile::getRepresentative()
{
	return &ms_RawIVFileRepresentative;
}

