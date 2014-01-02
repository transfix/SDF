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
/*****************************************************************************/// MRCFile.h: interface for the MRCFile class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CCV_MRCFILE_TYPE_H
#define CCV_MRCFILE_TYPE_H

#include "VolumeFileType.h"

typedef struct {
	
	int    nx;            // # of columns ( fastest changing in the map    
	int    ny;            // # of rows                                     
	
	int    nz;            // # of sections (slowest changing in the map    
	
	int    mode;          // data type
						  //      0 = image data in bytes
						  //      1 = image data in short integer
						  //      2 = image data in floats
						  //      3 = complex data in complex short integers
						  //      4 = complex data in complex reals          
	
	int    nxstart;       // number of first column in map (default = 0)   
	int    nystart;       // number of first row in map (default = 0)      
	int    nzstart;       // number of first ssection in map (default = 0) 
	
	int    mx;            // number of intervals along X                   
	int    my;            // number of intervals along Y                   
	int    mz;            // number of intervals along Z                   
	
	float  xlength;       // cell dimensions in X (angstrom)               
	float  ylength;       // cell dimensions in Y (angstrom)               
	float  zlength;       // cell dimensions in Z (angstrom)               
	
	float  alpha;         // cell angles between Y and Z                   
	float  beta;          // cell angles between X and Z                   
	float  gamma;         // cell angles between X and Y                   
	
	int    mapc;          // number of axis corresponding to columns (X)   
	int    mapr;          // number of axis corresponding to rows (Y)      
	int    maps;          // number of axis corresponding to sections (Z)  
	
	float  amin;          // minimum density value                         
	float  amax;          // maximum density value                         
	float  amean;         // mean density value                            
	
	int    ispg;          // space group number (0 for images)             
	int    nsymbt;        // # of bytes for symmetry operators             
	
	int    extra[29];     // user defined storage space                    
	
	float  xorigin;       // X phase origin                                
	float  yorigin;       // Y phase origin                                
	
	int    nlabl;         // # of labels being used in the MRC header      
	
	char   label[10][80]; // actual text labels                            
	
} OLDMrcHeader;

typedef struct {
	
	int    nx;            // # of columns ( fastest changing in the map    
	int    ny;            // # of rows                                     
	
	int    nz;            // # of sections (slowest changing in the map    
	
	int    mode;          // data type
						  //      0 = image data in bytes
						  //      1 = image data in short integer
						  //      2 = image data in floats
						  //      3 = complex data in complex short integers
						  //      4 = complex data in complex reals          
	
	int    nxstart;       // number of first column in map (default = 0)   
	int    nystart;       // number of first row in map (default = 0)      
	int    nzstart;       // number of first ssection in map (default = 0) 
	
	int    mx;            // number of intervals along X                   
	int    my;            // number of intervals along Y                   
	int    mz;            // number of intervals along Z                   
	
	float  xlength;       // cell dimensions in X (angstrom)               
	float  ylength;       // cell dimensions in Y (angstrom)               
	float  zlength;       // cell dimensions in Z (angstrom)               
	
	float  alpha;         // cell angles between Y and Z                   
	float  beta;          // cell angles between X and Z                   
	float  gamma;         // cell angles between X and Y                   
	
	int    mapc;          // number of axis corresponding to columns (X)   
	int    mapr;          // number of axis corresponding to rows (Y)      
	int    maps;          // number of axis corresponding to sections (Z)  
	
	float  amin;          // minimum density value                         
	float  amax;          // maximum density value                         
	float  amean;         // mean density value                            
	
	int    ispg;          // space group number (0 for images)             
	int    nsymbt;        // # of bytes for symmetry operators             
	
	int    extra[25];     // user defined storage space                    
	
	float  xorigin;       // X phase origin                                
	float  yorigin;       // Y phase origin                                
	float  zorigin;       // Z phase origin

	char   map[4];        // character string 'MAP '

	int    machst;        // machine stamp

	float  rms;           // rms deviation of map from mean density
	
	int    nlabl;         // # of labels being used in the MRC header      
	
	char   label[10][80]; // actual text labels                            
	
} MrcHeader;

class MRCFile : public VolumeFileType  
{
public:
	virtual ~MRCFile();

	virtual SimpleVolumeData* loadFile(const string& fileName);
	virtual bool checkType(const string& fileName);
	virtual bool saveFile(SimpleVolumeData* simpleVolumeData, const string& fileName, unsigned int variable=0);

	virtual string extension() { return "mrc"; };
	virtual string filter() { return "MRC files (*.mrc)"; };

	static MRCFile ms_MRCFileRepresentative;
	static VolumeFileType* getRepresentative();

protected:
	MRCFile();
	int finite(float f);

	bool readHeader(FILE* fp, int fileSize, bool* swap, SimpleVolumeData* simpleVolumeData);
	bool checkHeader(MrcHeader& header, unsigned int fileSize);
	void swapHeader(MrcHeader& header);
	bool interpretOldHeader(MrcHeader& header, unsigned int fileSize, SimpleVolumeData* simpleVolumeData);
	bool interpretNewHeader(MrcHeader& header, unsigned int fileSize, SimpleVolumeData* simpleVolumeData);
	bool fillHeader(MrcHeader& header, SimpleVolumeData* simpleVolumeData);

};

#endif
