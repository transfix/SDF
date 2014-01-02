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
/*     Author:              John Wiggins <prok@ices.utexas.edu>   2004-2005  */
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
/*****************************************************************************/
// MRCFile.cpp: implementation of the MRCFile class.
//
//////////////////////////////////////////////////////////////////////

#include "MRCFile.h"
#include "SimpleVolumeData.h"
#include "ByteSwapping.h"
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <cstring>

MRCFile MRCFile::ms_MRCFileRepresentative;

// XXX: This is UGLY. Windows does not have this function in its math library.
// This function is only called from the interpretXXXHeader functions
int MRCFile::finite(float f)
{
	if( f < -10e-10 || f > 10e10 )
		return 0;
	return 1;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
MRCFile::MRCFile()
{
	
}

MRCFile::~MRCFile()
{
	
}

/*
{
unsigned int dims[3]; dims[0] = dim1; dims[1] = dim2; dims[2] = dim3;
vol->setDimensions( dims );
}
vol->setNumberOfVariables(1);
vol->setData(0, m_SES);
vol->setType(0, SimpleVolumeData::FLOAT);
vol->setName(0, "skinRegion");
vol->setMinExtent(m_Min);
vol->setMaxExtent(m_Max);
*/
bool MRCFile::interpretNewHeader(MrcHeader& header, unsigned int fileSize, SimpleVolumeData* simpleVolumeData)
{
	unsigned int dims[3];
	dims[0] = header.nx;
	dims[1] = header.ny;
	dims[2] = header.nz;
	simpleVolumeData->setDimensions( dims );
	
	simpleVolumeData->setNumberOfVariables(1);
	
	simpleVolumeData->setName(0, "No Name");
	
	if( header.mode ==0 )
		simpleVolumeData->setType(0, SimpleVolumeData::UCHAR);
	else if( header.mode ==1 )
		simpleVolumeData->setType(0, SimpleVolumeData::USHORT);
	else if( header.mode ==2 )
		simpleVolumeData->setType(0, SimpleVolumeData::FLOAT);
	else
		return false;
	
	// make sure we aren't using garbage values
	float minExtent[3];
	if (!finite(header.xorigin) || !finite(header.yorigin) || !finite(header.zorigin)) {
		minExtent[0] = 0.0;
		minExtent[1] = 0.0;
		minExtent[2] = 0.0;
	}
	else {
		minExtent[0] = header.xorigin;
		minExtent[1] = header.yorigin;
		minExtent[2] = header.zorigin;
	}
	simpleVolumeData->setMinExtent(minExtent);
	
	float maxExtent[3];
	// we need to double check the meaning of xlength
	// (plus some extra paranoia)
	if (header.xlength<=0.0 || header.ylength<=0.0 || header.zlength<=0.0
		|| !finite(header.xlength) || !finite(header.ylength) || !finite(header.zlength)) {
		// hmm, this is wierd
		maxExtent[0] = minExtent[0] + dims[0]*1.0;
		maxExtent[1] = minExtent[1] + dims[1]*1.0;
		maxExtent[2] = minExtent[2] + dims[2]*1.0;
	}
	else {
		maxExtent[0] = minExtent[0] + header.xlength;
		maxExtent[1] = minExtent[1] + header.ylength;
		maxExtent[2] = minExtent[2] + header.zlength;
	}
	simpleVolumeData->setMaxExtent(maxExtent);
	
	// everything checks out, return true
	return true;
}

bool MRCFile::interpretOldHeader(MrcHeader& header, unsigned int fileSize, SimpleVolumeData* simpleVolumeData)
{
	unsigned int dims[3];
	dims[0] = header.nx;
	dims[1] = header.ny;
	dims[2] = header.nz;
	simpleVolumeData->setDimensions( dims );
	
	simpleVolumeData->setNumberOfVariables(1);
	
	simpleVolumeData->setName(0, "No Name");
	
	if( header.mode ==0 )
		simpleVolumeData->setType(0, SimpleVolumeData::UCHAR);
	else if( header.mode ==1 )
		simpleVolumeData->setType(0, SimpleVolumeData::USHORT);
	else if( header.mode ==2 )
		simpleVolumeData->setType(0, SimpleVolumeData::FLOAT);
	else
		return false;
	
	float minExtent[3];
	minExtent[0] = 0.0;
	minExtent[1] = 0.0;
	minExtent[2] = 0.0;
	simpleVolumeData->setMinExtent(minExtent);
	
	float maxExtent[3];
	// we need to double check the meaning of xlength
	// (plus some extra paranoia)
	if (header.xlength<=0.0 || header.ylength<=0.0 || header.zlength<=0.0
		|| !finite(header.xlength) || !finite(header.ylength) || !finite(header.zlength)) {
		// hmm, this is wierd
		maxExtent[0] = minExtent[0] + dims[0]*1.0;
		maxExtent[1] = minExtent[1] + dims[1]*1.0;
		maxExtent[2] = minExtent[2] + dims[2]*1.0;
	}
	else {
		maxExtent[0] = minExtent[0] + header.xlength;
		maxExtent[1] = minExtent[1] + header.ylength;
		maxExtent[2] = minExtent[2] + header.zlength;
	}
	simpleVolumeData->setMaxExtent(maxExtent);
	
	// everything checks out, return true
	return true;
}

bool MRCFile::checkHeader(MrcHeader& header, unsigned int fileSize)
{
	unsigned int sizes[] = {1, 2, 4};
	
	//qDebug("MrcFileImpl::checkHeader()");
	
	// check for the details we dont support
	if (header.mode<0 || header.mode>2) {
		// we dont support this type or MRC file for now		
		//qDebug("unsupported mrc file. (mode = %x)", header.mode);
		return false;
	}
	
	// check the fileSize
	if (sizes[header.mode]*header.nx*header.ny*header.nz + 1024 != fileSize) {
		// the size does not match the header information
		//qDebug("bad mrc file? (file size != size given in header)");
		return false;
	}
	
	// everything checks out, return true
	return true;
}

void MRCFile::swapHeader(MrcHeader& header)
{
	//qDebug("swap header");
	swapByteOrder(header.nx);
	swapByteOrder(header.ny);
	
	swapByteOrder(header.nz);
	
	//qDebug("header.mode = %x", header.mode);
	swapByteOrder(header.mode);
	//qDebug("header.mode = %x", header.mode);
	
	swapByteOrder(header.nxstart);
	swapByteOrder(header.nystart);
	swapByteOrder(header.nzstart);
	
	swapByteOrder(header.mx);
	swapByteOrder(header.my);
	swapByteOrder(header.mz);
	
	swapByteOrder(header.xlength);
	swapByteOrder(header.ylength);
	swapByteOrder(header.zlength);
	
	swapByteOrder(header.alpha);
	swapByteOrder(header.beta);
	swapByteOrder(header.gamma);
	
	swapByteOrder(header.mapc);
	swapByteOrder(header.mapr);
	swapByteOrder(header.maps);
	
	swapByteOrder(header.amin);
	swapByteOrder(header.amax);
	swapByteOrder(header.amean);
	
	swapByteOrder(header.ispg);
	swapByteOrder(header.nsymbt);
	
	swapByteOrder(header.extra, 25);
	
	swapByteOrder(header.xorigin);
	swapByteOrder(header.yorigin);
	swapByteOrder(header.zorigin);
	
	swapByteOrder(header.rms);
	
	swapByteOrder(header.nlabl);
	
}

bool MRCFile::readHeader(FILE* fp, int fileSize, bool* swap, SimpleVolumeData* simpleVolumeData)
{
	MrcHeader header;
	fread( &header, sizeof(MrcHeader), 1, fp );
	
	if (!(header.map[0]=='M' && header.map[1]=='A' && header.map[2]=='P')) {
		// not the new MRC style, we must try to guess the
		// endianness
		
		// first try not swapping
		if (checkHeader(header, fileSize)) {
			*swap = false;
			return interpretOldHeader(header, fileSize, simpleVolumeData);
		}
		else {
			// swap and try again
			swapHeader(header);
			if (checkHeader(header, fileSize)) {
				*swap = true;
				return interpretOldHeader(header, fileSize, simpleVolumeData);
			}
			else {
				// we dont support wierd or exotic endianness
				return false;
			}
		}
	}
	// Nobody seems to agree about the meaning of the machine stamp,
	// so swap again if the header doesn't check out.
	if (!checkHeader(header, fileSize)) {
		// we change our mind about swapping
		//m_MustSwap = !m_MustSwap;
		*swap = true;
		// and swap the header again
		swapHeader(header);
	}
	if (!checkHeader(header, fileSize)) {
		return false;
	}
	
	
	return interpretNewHeader(header, fileSize, simpleVolumeData);
}

SimpleVolumeData* MRCFile::loadFile(const string& fileName)
{
	SimpleVolumeData* simpleVolumeData = new SimpleVolumeData(128, 128, 128);
	
	int fileSize = 0;
	{
		FILE *fp=fopen(fileName.c_str(), "rb");
		if (!fp) {
			printf("Error: could not open file\n");
			delete simpleVolumeData; simpleVolumeData = 0;
			return 0;
		}
		int start = ftell(fp );
		fseek(fp, 0, SEEK_END);
		int end = ftell( fp );
		fileSize = end-start;
		fclose(fp );
		if( start<0 || end <0 || fileSize<1 ) 
		{
			printf("Error: in determining file size\n");
			delete simpleVolumeData; simpleVolumeData = 0;
			return 0;
		}
	}
	FILE *fp=fopen(fileName.c_str(), "rb");
	
	// check to make sure the file exists
	if (!fp) {
		printf("Error: could not open file\n");
		delete simpleVolumeData; simpleVolumeData = 0;
		return 0;
	}
	
	// read the header
	// 
	bool swap = false;
	if( !readHeader(fp, fileSize, &swap, simpleVolumeData) ) 
	{
		delete simpleVolumeData; simpleVolumeData = 0;
		return 0;
	}
	
	int numVerts =  simpleVolumeData->getWidth()*
		simpleVolumeData->getHeight()*
		simpleVolumeData->getDepth();
	int dataSize =	simpleVolumeData->getTypeSize(0)*
		numVerts;
	
	if( dataSize < 1 )
	{
		delete simpleVolumeData; simpleVolumeData = 0;
		return 0;
	}
	unsigned char* data = new unsigned char[dataSize];
	fread(data, dataSize, 1, fp);
	if (swap) {
		switch(simpleVolumeData->getTypeSize(0)) {
		case 1:
			break;
		case 2:
			swapByteOrder((unsigned short *)data, numVerts);
			break;
		case 4:
			swapByteOrder((float *)data, numVerts);
			break;
		case 8:
			swapByteOrder((double *)data, numVerts);
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

bool MRCFile::checkType(const string& fileName)
{
	return false;
}

bool MRCFile::fillHeader(MrcHeader& header, SimpleVolumeData* simpleVolumeData)
{
	// fill in the header's fields
	header.nx = simpleVolumeData->getWidth();
	header.ny = simpleVolumeData->getHeight();
	header.nz = simpleVolumeData->getDepth();
	
	switch (simpleVolumeData->getType(0))
	{
	case SimpleVolumeData::UCHAR:
		header.mode = 0;
		break;
	case SimpleVolumeData::USHORT:
		header.mode = 1;
		break;
	case SimpleVolumeData::FLOAT:
		header.mode = 2;
		break;
	default:
		// wha???
		// we shouldn't get here. (famous last words, I know)
		return false;
		break;
	}
	
	// start coord, defaults to (0,0,0)
	header.nxstart = 0;   
	header.nystart = 0;      
	header.nzstart = 0; 
	
	// the dimensions again
	header.mx = simpleVolumeData->getWidth();
	header.my = simpleVolumeData->getHeight();
	header.mz = simpleVolumeData->getDepth();
	
	// dimensions of a cell (span) 
	// (supposed to be in angstroms, but no guarantees)
	header.xlength = simpleVolumeData->getMaxX() - simpleVolumeData->getMinX();
	header.ylength = simpleVolumeData->getMaxY() - simpleVolumeData->getMinY();
	header.zlength = simpleVolumeData->getMaxZ() - simpleVolumeData->getMinZ();
	
	// cell angles, all 90 deg
	header.alpha = 90.0;
	header.beta = 90.0;
	header.gamma = 90.0;
	
	// axis order
	header.mapc = 1; // number of axis corresponding to columns (X)
	header.mapr = 2; // number of axis corresponding to rows (Y)
	header.maps = 3; // number of axis corresponding to sections (Z)
	
	// min, max and mean... just put 0.0
	header.amin = 0.0; // minimum density value
	header.amax = 0.0; // maximum density value
	header.amean = 0.0; // mean density value
	
	header.ispg = 0; // space group number (0 for images)
	header.nsymbt = 0; // # of bytes for symmetry operators
	
	memset(header.extra, 0, 25*sizeof(int)); // user defined storage space
	
	// mesh origin
	header.xorigin = simpleVolumeData->getMinX(); // X phase origin
	header.yorigin = simpleVolumeData->getMinY(); // Y phase origin
	header.zorigin = simpleVolumeData->getMinZ(); // Z phase origin
	
	// character string 'MAP '
	header.map[0] = 'M';
	header.map[1] = 'A';
	header.map[2] = 'P';
	header.map[3] = ' ';
	
	// machine stamp
	if (isLittleEndian()) {
		header.machst = 0x44410000;
		// swap it to big endian
		swapByteOrder(header.machst);
	}
	else
		header.machst = 0x11110000;
	
	header.rms = 0.0; // rms deviation of map from mean density
	header.nlabl = 1; // # of labels being used in the MRC header
	
	// zero the labels
	for (int i=0; i < 10; i++)
		memset(header.label[i], 0, 80);
	// fill in the first label
	strcpy(header.label[0], "Created by TexMol");
	
	return true;
}

bool MRCFile::saveFile(SimpleVolumeData* simpleVolumeData, const string& fileName, unsigned int variable)
{
	FILE *fp=fopen(fileName.c_str(), "wb");
	
	// check to make sure the file exists
	if (!fp) {
		printf("Error: could not open file\n");
		delete simpleVolumeData; simpleVolumeData = 0;
		return false;
	}

	MrcHeader header;
	
	// fill in the header
	if( !fillHeader(header, simpleVolumeData) ) 
	{
		fclose( fp );
		return false;
	}
	
	//if (isBigEndian())
	//	swapHeader(header);
	
	fwrite( &header, sizeof( MrcHeader ), 1, fp );

	int dataSize = simpleVolumeData->getDataSize(0);
	if( dataSize < 1 )
	{
		fclose( fp );
		return false;
	}
	fwrite( simpleVolumeData->getData(0), dataSize, 1, fp );

	fclose( fp );
	return true;
}

VolumeFileType* MRCFile::getRepresentative()
{
	return &ms_MRCFileRepresentative;
}
