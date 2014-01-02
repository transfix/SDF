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
/*****************************************************************************/// DXFile.h: interface for the DXFile class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CCV_DX_FILE_H
#define CCV_DX_FILE_H

#include "VolumeFileType.h"

class DXFile : public VolumeFileType  
{
public:
	virtual ~DXFile();

	virtual SimpleVolumeData* loadFile(const string& fileName);
	virtual bool checkType(const string& fileName);
	virtual bool saveFile(SimpleVolumeData* simpleVolumeData, const string& fileName, unsigned int variable=0);

	virtual string extension() { return "dx"; };
	virtual string filter() { return "DX files (*.dx)"; };

	static DXFile ms_DXFileRepresentative;
	static VolumeFileType* getRepresentative();

protected:
	DXFile();
	bool isCommentOrEmpty( const char* line );
	void tryToGetDimensions( bool* dimsSet, unsigned int* dims, const char* line );
	void tryToGetOrigin( bool* originSet, float* origin, const char* line );
	void tryToGetSpans( int* deltasFound, float* span, const char* line );
	void tryToGetDataHeader(bool* dataHeaderFound, unsigned int* dims, const char* line);

};

#endif 
