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
/*****************************************************************************/

#include "Server.h"
#include "SurfaceDataManager/SurfaceData.h"
#include "BallAndStickDataManager/BallAndStickDataManager.h"
#include "DockingManager.h"
#include "SimpleVolumeData.h"
#include "moleculevizmainwindow.h"
#include "ColorTable.h"
#include "GroupOfAtoms.h"
#include "Atom.h"
#include "parserPDBtoGOA.h"
#include "GOAFileIO.h"
#include "MultiContour.h"
#include "VolumeDataManager/VolumeData.h"
#include "GeometryLoader.h"
#include "VolumeLoader.h"
#include "GOALoader.h"
#include "SimpleVolumeDataIsocontourer.h"
#include "UsefulMath/LinearAlgebra.h"
#include "BlurMapsDataManager.h"
#include "SurfaceAtomExtractor.h"
#include "SkinRegion2.h"
#include "contour.h"
#include "DownloadPDB.h"
#include "AreaVolume.h"
#include "sdfLib.h"
#include "MolecularCharacteristics.h"
#include "MergeVolumes.h"
#include "MoleculeMorph.h"

#include "DirectToGridSummationModule.h"
#include "GaussianKernel.h"
#include "UniformOutputGrid.h"

#include "Matrix.h"
#include "Vector.h"
#include "Quaternion.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <qfile.h>
#include <qtextstream.h>

extern void getContourSpectrum(unsigned char* uchar_data, int type, int* dim, int array_size, float* isoval , float* area, float* min_vol, float* max_vol, float* gradient, float* span);

Server::Server()
{
}

Server::~Server()
{

}
bool Server::depthColor( int argc, char* argv[] )
{
	// Take as input a rawiv file and create a depth colored rawv.
	// It needs a color map also. We are currently using a default map

	if( argc != 7 )
	{
		printUsage();
		return false;
	}

	char inputFileName[256];
	char outputFileName[256];
	int dim1, dim2, dim3;
	
	strcpy( inputFileName, argv[2] );
	strcpy( outputFileName, argv[3] );
	dim1 = atoi(argv[4]); dim2 = atoi(argv[5]); dim3 = atoi(argv[6]);

	int colorMapSize = 256; // new rover supports only 256  - be careful SKVINAY !
	double dColorMap[256*4];
	
	int i, j;
	for( i=0; i<256; i++ )
	{
		dColorMap[i*4+0] = 1;
		dColorMap[i*4+1] = 0;
		dColorMap[i*4+2] = 0;
		dColorMap[i*4+3] = 1;
	}

	double r=1.0, g=1.0, b=1.0;
	double oldr = 1.0, oldg = 1.0, oldb = 1.0;
	for( i=0; i<8; i++ )
	{
		for( j=0; j<32; j++ )
		{
			dColorMap[(i*32 + j)*4 +0] = (r*j/32.0 + oldr*(32-j)/32);
			dColorMap[(i*32 + j)*4 +1] = (g*j/32.0 + oldg*(32-j)/32);
			dColorMap[(i*32 + j)*4 +2] = (b*j/32.0 + oldb*(32-j)/32);
			dColorMap[(i*32 + j)*4 +3] = 1;
		}
		printf("%d %f %f %f\n", i, r, g, b );
		oldr = r; oldg = g; oldb = b;
		r = rand()/((double)(RAND_MAX)+1); b = rand()/((double)(RAND_MAX)+1); g = rand()/((double)(RAND_MAX)+1);
	}
	SimpleVolumeData* simpleVolumeFile = 0;
	{
		VolumeLoader* volumeLoader = new VolumeLoader();
		simpleVolumeFile = volumeLoader->loadFile( inputFileName );
		delete volumeLoader;
		if( !simpleVolumeFile ) return false;
	}

	//////////////  create depth colored volume ////////////////////
	SimpleVolumeData* sData = 0;
	sData = simpleVolumeFile->createDepthColoredVolume(dColorMap, colorMapSize-1 );
	if( !sData ) 
	{
		delete simpleVolumeFile;
		return false;
	}
	////////////////////////////////////////////////////////////////

	/////////////// write to file if needed ////////////////////////
	bool ret = false;
	{
		VolumeLoader* volumeLoader = new VolumeLoader();
		ret = volumeLoader->saveFile( outputFileName, sData );
		delete volumeLoader;
	}
	/////////////////////////////////////////////////////////////////

	delete simpleVolumeFile;
	delete sData;
	return ret;
}

bool Server::blur( int argc, char* argv[] )
{
	// Create a blur map of an input pdb or pqr file. Should we also support pts files ?
	// Supports both electron density and hydrophobicity
	// Input parameters:
	//		char* input filename
	//		char* output filename
	//		int dim1, dim2, dim3
	//		double blobby factor
	//		int density type ( 0 - electron density, 1 - hydrophobicity )
	//		char* colormap
	//		int color type ( 0 - , 1 - , 2 - )
	//		int gap - not very nice to use as it destroys meaning of origin, span etc. Only use sparingly for viz
	//	Side effects:
	//		creates RAWIV or RAV5 volume 
	//	Return value:
	//		bool indicates success or failure of function call
	//	usage:
	//			bool getVolume( const char* pdbOrPqrFileName, const char* volFileName, 
	//				 int dim1, int dim2, int dim3, int densityType,
	//				 bool writeRawV, double blob, unsigned int color, const char* cmapFile, int gap, int radiusType, unsigned int level );
	//  radiusType == 0 -> vdv radius, == 1 -> solvent enlarged radius
	//
	// the 16 to 27 parameters can set to transform the molecule: (cx, cy, cz, nx, ny, nz) old, (cx, cy, cz, nx, ny, nz) new
	if( argc != 15 && argc != 27) 
	{
		printUsage();
		return false;
	}

	char inputFileName[256];
	char outputFileName[256];
	int dim1, dim2, dim3;
	PDBParser::GroupOfAtoms::FUNCTIONS densityType;
	bool writeRawV;
	double blob;
	PDBParser::GroupOfAtoms::GOA_TYPE colorLevel;
	char cmapFile[256];
	int gap;
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	int level;

	strcpy( inputFileName, argv[2] );
	strcpy( outputFileName, argv[3] );
	dim1 = atoi(argv[4]); dim2 = atoi(argv[5]); dim3 = atoi(argv[6]);
	if( !PDBParser::GroupOfAtoms::intToFunctionType( &densityType, atoi(argv[7]) ) ) return false;
	if( strcmp(argv[8], "true") == 0 ) 
		writeRawV = true;
	else if( strcmp(argv[8], "false") == 0 ) 
		writeRawV = false;
	else
	{
		printUsage();
		return false;
	}
	blob = atof(argv[9]);
	if( !PDBParser::GroupOfAtoms::intToGOAType( &colorLevel, atoi(argv[10]) ) ) return false;
	strcpy( cmapFile, argv[11] );
	gap = atoi(argv[12]); 
	if( !PDBParser::GroupOfAtoms::intToRadiusType( &radiusType, atoi(argv[13]) ) ) return false;
	level = atoi(argv[14] );

	SimpleVolumeData* sData = 0;
	///////////// transform molecule if required////////
	if( argc == 27 )
	{
		CCVOpenGLMath::Matrix transformation;
		{
			CCVOpenGLMath::Vector old_center(atof(argv[15]), atof(argv[16]), atof(argv[17]), 1);
			CCVOpenGLMath::Vector old_normal(atof(argv[18]), atof(argv[19]), atof(argv[20]), 0);

			CCVOpenGLMath::Vector new_center(atof(argv[21]), atof(argv[22]), atof(argv[23]), 1);
			CCVOpenGLMath::Vector new_normal(atof(argv[24]), atof(argv[25]), atof(argv[26]), 0);

			// shift center to new one.
			// rotate to align quaternions
			// for each mol, we need a predefined center and normal.
			// translate mol_center to input_center
			// now to rotate:
			//     cross product gives axis of rotation
			//     dot product gives the angle
			old_normal.normalize();
			new_normal.normalize();
			CCVOpenGLMath::Vector axis_of_rotation = old_normal.cross( new_normal );
			axis_of_rotation.normalize();
			double angle_of_rotation = acos(new_normal.dot( old_normal ));

			float w = cos( angle_of_rotation / 2.0 );
			float x = axis_of_rotation[0] * sin( angle_of_rotation / 2.0 );
			float y = axis_of_rotation[1] * sin( angle_of_rotation / 2.0 );
			float z = axis_of_rotation[2] * sin( angle_of_rotation / 2.0 );
			CCVOpenGLMath::Quaternion quaternion(w, x, y, z);
			CCVOpenGLMath::Matrix rotation = quaternion.buildMatrix();

			CCVOpenGLMath::Vector tr_old_center = rotation*old_center;
			CCVOpenGLMath::Matrix translation = CCVOpenGLMath::Matrix::translation(new_center - tr_old_center);

			transformation = rotation.preMultiplication(translation);
			//transformation = rotation;
		}

		sData = BlurMapsDataManager::getVolume( inputFileName, outputFileName, 
			    dim1, dim2, dim3, densityType,
				writeRawV, blob, colorLevel, cmapFile, gap, radiusType, level, &transformation);

	}
	////////////////////////////////////////////////////
	else
	{
		sData = BlurMapsDataManager::getVolume( inputFileName, outputFileName, 
			    dim1, dim2, dim3, densityType,
				writeRawV, blob, colorLevel, cmapFile, gap, radiusType, level, 0);
	}

	bool ret;
	if( sData ) 
		ret = true;
	else
		ret = false;
	delete sData;
	return ret;
/*{
		GOALoader* gLoader = new GOALoader();
		PDBParser::GroupOfAtoms* molecule = gLoader->loadFile( argv[2] );
		delete gLoader;
		if( !molecule ) return false;

		int n = 0;
		molecule->getNumberOfAtomsRecursive(&n);
		if( n <= 0 ) return false;

		double blobbiness = atof(argv[9]);
		if( blobbiness >= 0 ) return false;

		double error = 1e-10;
		
		PDBParser::GroupOfAtoms::FUNCTIONS function = PDBParser::GroupOfAtoms::ELECTRON_DENSITY;
		if( !PDBParser::GroupOfAtoms::intToFunctionType( &function, atoi(argv[7]) ) ) return false;

		PDBParser::GroupOfAtoms::RADIUS_TYPE radius_type = PDBParser::GroupOfAtoms::VDW_RADIUS;
		if( !PDBParser::GroupOfAtoms::intToRadiusType( &radius_type, atoi(argv[13]) ) ) return false;

		int GOA_level = atoi(argv[14] );
		
		double min[3]; min[0] = min[1] = min[2] = 1e20;
		double max[3]; max[0] = max[1] = max[2] = 1e-20;
		double maxRadius = -1;

		double* points = new double[n*3];
		double* radii = new double[n];
		double* weights = new double[n];

		molecule->getAttributes( points, radii, min, max, weights, &maxRadius, function, GOA_level, radius_type, n );
		CCVSummationModule::Kernel* kernel = new CCVSummationModule::GaussianKernel( blobbiness, error );

		unsigned int dimensions[3];
		dimensions[0] = atoi(argv[4]); dimensions[1] = atoi(argv[5]); dimensions[2] = atoi(argv[6]);
		float* output = new float[dimensions[0]*dimensions[1]*dimensions[2]];
		float origin[3]; origin[0] = min[0]; origin[1] = min[1]; origin[2] = min[2];
		float span[3];
		span[0] = (max[0] - min[0]) / ((double)(dimensions[0]-1));
		span[1] = (max[1] - min[1]) / ((double)(dimensions[1]-1));
		span[2] = (max[2] - min[2]) / ((double)(dimensions[2]-1));

		CCVSummationModule::OutputGrid *  outputGrid = new CCVSummationModule::UniformOutputGrid( output, origin, span, dimensions );

		CCVSummationModule::SummationModule* summationModule = new CCVSummationModule::DirectToGridSummationModule( points, radii, weights, n, kernel, outputGrid );
		if( !summationModule->sum() ) 
		{
			if( points ) { delete []points; points = 0; }
			if( radii ) { delete []radii; radii = 0; }
			if( weights ) { delete []weights; weights = 0; }
			if( output ) { delete []output; output = 0; }
			return false;
		}
		{
			SimpleVolumeData* vol = new SimpleVolumeData(dimensions);
			vol->setDimensions( dimensions );
			vol->setNumberOfVariables(1);
			vol->setData(0, output);
			vol->setType(0, SimpleVolumeData::FLOAT);
			vol->setName(0, "TexMols blurring");
			float fmin[3], fmax[3];
			fmin[0] = min[0]; fmin[1] = min[1]; fmin[2] = min[2]; 
			fmax[0] = max[0]; fmax[1] = max[1]; fmax[2] = max[2]; 
			vol->setMinExtent(fmin);
			vol->setMaxExtent(fmax);

			VolumeLoader* volumeLoader = new VolumeLoader();
			volumeLoader->saveFile("opt.rawiv", vol );
			delete volumeLoader;

			delete vol;
		}

		// simpleVolumeData should delete the 'output' data
		if( points ) { delete []points; points = 0; }
		if( radii ) { delete []radii; radii = 0; }
		if( weights ) { delete []weights; weights = 0; }
	}
	return true;
	*/
}

void Server::printUsage()
{
	printf("Usage: any one of the following, all corresponding paramters are required. Consult manual for more help\n");
	printf("1.  MoleculeViz\n");
	printf("2.  MoleculeViz [-blur <const char* pdbOrPqrFileName> <const char* volFileName> \n<int dim1> <int dim2> <int dim3> <int densityType> \n<bool writeRawV> <double blob> <unsigned int color> <const char* cmapFile> <int gap> <int radiusType> <unsigned int level>]\n");
	printf("WRONG , CHANGE 3.  MoleculeViz [-setcurvature <int createIsosurface> <const char* inputpdbOrPqrFileName> <int dim1> <int dim2> <int dim3> <double blob> <const char* inputRawSurfaceFileName> <const char* outputMeanRawSurfaceFileName> <const char* outputGaussianRawSurfaceFileName> <const char* outputCurvatureRawSurfaceFileName> <int numberOfGridDivisions> <double maxFunctionError> <int radiusType>]\n");
	printf("4.  MoleculeViz [-outGridPositions <const char* pdbOrPqrFileName> <unsigned int N1> <unsigned int N2> <unsigned int N3> <double extraSpace>]\n");
	printf("5.  MoleculeViz [-classifyPoints <const char* inputPtsFileName> <const char* outputPtsFileName> ]\n");
	printf("6.  MoleculeViz [-growOut <const char* inputRawFileName> <const char* outputRawFileName> <double extent>]\n");
	printf("7.  MoleculeViz [-getSurfaceFromVolume <const char* volumeFileName> <const char* outputSurfaceFileName> <double isovalue>]\n");
	printf("8.  MoleculeViz [-getSurfaceFromPDB <const char* pdbOrPqrFileName> <const char* outputSurfaceFileName> <double isovalue> <int dim1> <int dim2> <int dim3> <double blobbiness>]\n");
	printf("9.  MoleculeViz [-evolve <const char* pdbOrPqrFileName> <const char* changesFile> <const char* pdbOrPqrFileName>]\n");
	printf("10. MoleculeViz [-writePDB <const char* inputFileName> <const char* outputFileName> <int outputLevel>]\n");
	printf("11. MoleculeViz [-writeGOA <const char* inputFileName> <const char* outputFileName>]\n");
	printf("12. MoleculeViz [-depthColor <const char* inputFileName> <const char* outputFileName> <int dim1> <int dim2> <int dim3>]\n");
	printf("13. MoleculeViz [-getHydrophobicityOnSurface <const char* inputPDBFileName> <const char* inputSurfaceFileName> <int dim1> <int dim2> <int dim3> <double blobbiness> <const char* outputValuesFileName>]\n");
	printf("14. MoleculeViz [-createSkinRegion <const char* inputPDBFileName> <const char* inputSurfaceFileName> <int dim1> <int dim2> <int dim3> <double probeRadius> <int radiusType>]\n");
	printf("15. MoleculeViz [-correlate <const char* hydroFile> <const char* meanCurvFile> <const char* gausCurvFile> <double probeRadius> <bool discretizeHydro> <bool discretizeCurv> <double maxDiscreteVal> <double minDiscreteVal> <int rangeToCorrelate = {0->all,<0,>0} >]\n");
	printf("16. MoleculeViz [-getArea <const char* surfaceFileName>]\n");
	printf("17. MoleculeViz [-getPatchAreas <const char* surfaceFileName> <const char* functionValuesFileName> <double isovalue> <const char* outputFileName>]\n");
	printf("18. MoleculeViz [-populateSAS <const char* surfaceFileName> <const char* functionValuesFileName> <double isovalue> <const char* outputFileName>]\n");
	printf("19. MoleculeViz [-addVolumes <const char* volume1FileName> <const char* volume2FileName> <const char* volume3FileName> <double scale> <double sum>]\n");
	printf("20. MoleculeViz [-convert <const char* inputFileName> <const char* outputFileName>]\n");
	printf("21. MoleculeViz [-getSurfaceAtoms <const char* inputFileName> <const char* outputFileName>]\n");
	printf("22. MoleculeViz [-getContourStats <const char* inputFileName> <const char* outputFileName>]\n");
	printf("23. MoleculeViz [-downloadPDB <const char* pdbID> <const char* outputFileName>]\n");
	printf("24. MoleculeViz [-getSignDistanceFunction <const char* inputFileName> <const char* outputFileName> <int size> <int flipNormals>]\n");
	printf("25. MoleculeViz [-getMolecularCharacteristics <const char* inputGOAFileName> <const char* elecFileName> <const char* outputFileName> <int numberOfSpheres> opt: <float* active site> <double radius1> .. <double radius_n_of_spheres>]\n" );
	printf("26. MoleculeViz [-expandMolecule <const char* inputGOAFileName> <const char* outputGOAFileName>]\n" );
	printf("27. MoleculeViz [-mergeGeometry <const char* inputSurfaceFileName1> <const char* inputSurfaceFileName2> <const char* outputSurfaceFileName>]\n" );
	printf("28. MoleculeViz [-mergeVolumes <const char* inputSurfaceFileName1> <const char* inputSurfaceFileName2> <const char* outputSurfaceFileName> <optional: original and new center and normals>]\n" );
	printf("29. MoleculeViz [-writeTorsionAngles <const char* inputGOAFileName> <const char* outputTorsionAnglesFileName>]\n" );
	printf("30. MoleculeViz [-morph <const char* inputGOAFileName1> <const char* inputGOAFileName2> <const char* outputGOAFileNamePrefix> <double resolution> <int maxSteps>]\n" );
	printf("31. MoleculeViz [-getElecOnSurface <const char* inputElecVolume> <const char* inputSurfaceFileName>  <const char* outputPosSurfaceFileName>  <const char* outputNegSurfaceFileName>  <const char* outputNeuSurfaceFileName> <double negCutoff> <double posCutoff>]\n" );
	printf("32. MoleculeViz [-getCurvaturesOnSurface <const char* inputSurfaceFileName> <const char* outputFileName>]\n");
}

bool Server::execute( int argc, char* argv[], MoleculeVizMainWindow* mWindow )
{
	if( argc < 2 ) return false;

	if( strcmp(argv[1],"-help") == 0 )
	{
		printUsage();
		return true;
	}	
	if( strcmp(argv[1],"-h") == 0 )
	{
		printUsage();
		return true;
	}
	if( strcmp(argv[1],"-blur") == 0 )
		return blur(argc, argv);
	if( strcmp(argv[1],"-setcurvature") == 0 )
		return setCurvature(argc, argv);
	if( strcmp(argv[1],"-outGridPositions") == 0 )
		return outGridPositions(argc, argv);
	if( strcmp(argv[1],"-classifyPoints") == 0 )
		return classifyPoints(argc, argv);
	if( strcmp(argv[1],"-growOut") == 0 )
		return growOut(argc, argv);
	if( strcmp(argv[1],"-getSurfaceFromVolume") == 0 )
		return getSurfaceFromVolume(argc, argv);
	if( strcmp(argv[1],"-getSurfaceFromPDB") == 0 )
		return getSurfaceFromPDB(argc, argv);
	if( strcmp(argv[1],"-evolve") == 0 )
		return evolve(argc, argv);
	if( strcmp(argv[1],"-writePDB") == 0 )
		return writePDB(argc, argv);
	if( strcmp(argv[1],"-writeGOA") == 0 )
		return writeGOA(argc, argv);
	if( strcmp(argv[1],"-depthColor") == 0 )
		return depthColor(argc, argv);
	if( strcmp(argv[1],"-getMaxDistanceFromPoint") == 0 )
		return getMaxDistanceFromPoint(argc, argv);
	if( strcmp(argv[1],"-getHydrophobicityOnSurface") == 0 )
		return getHydrophobicityOnSurface(argc, argv);
	if( strcmp(argv[1],"-createSkinRegion") == 0 )
		return createSkinRegion(argc, argv);
	if( strcmp(argv[1],"-correlate") == 0 )
		return correlate(argc, argv);
	if( strcmp(argv[1],"-getArea") == 0 )
		return getArea(argc, argv);
	if( strcmp(argv[1],"-getPatchAreas") == 0 )
		return getPatchAreas(argc, argv);
	if( strcmp(argv[1],"-populateSAS") == 0 )
		return populateSAS(argc, argv);
	if( strcmp(argv[1],"-addVolumes") == 0 )
		return addVolumes(argc, argv);
	if( strcmp(argv[1],"-convert") == 0 )
		return convert(argc, argv);
	if( strcmp(argv[1],"-getSurfaceAtoms") == 0 )
		return getSurfaceAtoms(argc, argv);
	if( strcmp(argv[1],"-getContourStats") == 0 )
		return getContourStats(argc, argv);
	if( strcmp(argv[1],"-getSignDistanceFunction") == 0 )
		return getSignDistanceFunction(argc, argv);
	if( strcmp(argv[1],"-getMolecularCharacteristics") == 0 )
		return getMolecularCharacteristics(argc, argv);
	if( strcmp(argv[1],"-expandMolecule") == 0 )
		return expandMolecule(argc, argv);
	if( strcmp(argv[1],"-mergeGeometry") == 0 )
		return mergeGeometry(argc, argv);	
	if( strcmp(argv[1],"-mergeVolumes") == 0 )
		return mergeVolumes(argc, argv);	
	if( strcmp(argv[1],"-writeTorsionAngles") == 0 )
		return writeTorsionAngles(argc, argv);	
	if( strcmp(argv[1],"-morph") == 0 )
		return morph(argc, argv);
	if( strcmp(argv[1],"-getElecOnSurface") == 0 )
		return getElecOnSurface(argc, argv);
	if( strcmp(argv[1],"-getCurvaturesOnSurface") == 0 )
		return getCurvaturesOnSurface(argc, argv);

	if( strcmp(argv[1],"-addData") == 0 )
		return addNewDataSet( argc, argv, mWindow );
	if( strcmp(argv[1],"-splitView") == 0 )
		return splitView( argc, argv, mWindow );
	if( strcmp(argv[1],"-deleteData") == 0 )
		return deleteData( argc, argv, mWindow );
	if( strcmp(argv[1],"-deletePrevData") == 0 )
		return deletePrevData( argc, argv, mWindow );
	if( strcmp(argv[1],"-deleteAllData") == 0 )
		return deleteAllData( argc, argv, mWindow );
	if( strcmp(argv[1],"-setVisible") == 0 )
		return setVisible( argc, argv, mWindow );
	if( strcmp(argv[1],"-setVisiblePrev") == 0 )
		return setVisiblePrev( argc, argv, mWindow );
	if( strcmp(argv[1],"-saveImage") == 0 )
		return saveImage( argc, argv, mWindow );
	if( strcmp(argv[1],"-setGridVisible") == 0 )
		return setGridVisible( argc, argv, mWindow );
	if( strcmp(argv[1],"-downloadPDB") == 0 )
		return downloadPDB(argc, argv, mWindow);

	return false;
}

bool Server::outGridPositions( int argc, char* argv[] )
{
	if( argc != 7 )
	{
		printUsage();
		return false;
	}

	char pdbOrPqrFileName[256];
	unsigned int N1;
	unsigned int N2;
	unsigned int N3;
	double extraSpace;

	strcpy( pdbOrPqrFileName, argv[2] );
	N1 = atoi(argv[3]); N2 = atoi(argv[4]); N3 = atoi(argv[5]);
	extraSpace = atof( argv[6] );
	DockingManager::outGridPositionsInt( pdbOrPqrFileName, NULL, N1, N2, N3, extraSpace );
	//DockingManager::outGridPositions( pdbOrPqrFileName, N1, N2, N3, extraSpace );

	return true;
}

bool Server::getSurfaceFromVolume( int argc, char* argv[] )
{
	if( argc != 5 )
	{
		printUsage();
		return false;
	}

	char volumeFileName[256];
	char surfaceFileName[256];
	double isovalue = 1;

	strcpy( volumeFileName, argv[2] );
	strcpy( surfaceFileName, argv[3] );
	isovalue = atof( argv[4] );


	///////// read in the volume ////////////////////
	SimpleVolumeData* sData = 0;
	{
		VolumeLoader* vLoader = new VolumeLoader();
		sData = vLoader->loadFile( volumeFileName );
		delete vLoader;
		if( !sData )return false;
	}
	/////////////////////////////////////////////////




	Geometry* geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);

	if( geometry == 0 ) 
	{
		delete sData;
		return false;
	}
	///////////////////////////////////////////////////



	///////// save the isocontour /////////////////////
	GeometryLoader* geometryLoader = new GeometryLoader();
	if( !geometryLoader->saveFile( surfaceFileName, "Rawnc files (*.rawnc)", geometry) )
	{
		delete geometry;
		delete sData;
		return false;
	}
	delete geometryLoader;
	///////////////////////////////////////////////////

	delete geometry;
	delete sData;
	return true;
}

bool Server::evolve( int argc, char* argv[] )
{
	if( argc != 5 )
	{
		printUsage();
		return false;
	}

	char pdbOrPqrFileName[256];
	char changesFile[256];
	char outputPdbOrPqrFileName[256];

	strcpy( pdbOrPqrFileName, argv[2] );
	strcpy( changesFile, argv[3] );
	strcpy( outputPdbOrPqrFileName, argv[4] );

	BallAndStickDataManager b;
	b.evolve( pdbOrPqrFileName, changesFile, outputPdbOrPqrFileName );

	return true;
}

bool Server::setCurvature( int argc, char* argv[] )
{
	if( argc < 3 )
	{
		printUsage();
		return false;
	}

	char inputPQRorPDBFileName[256];
	char inputRawSurfaceFileName[256];
	char outputMeanRawSurfaceFileName[256];
	char outputGaussianRawSurfaceFileName[256];
	char outputVolumeFileName[256];
	char curvatureFileName[256];
	int dim1, dim2, dim3;
	double blob, isovalue;
	int createIsosurface;
	int numberOfGridDivisions;
	double maxFunctionError;
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	int level;

	createIsosurface = atoi(argv[2]);
	if( createIsosurface == 0 )
	{
		if( argc != 16 )
		{
			printUsage();
			return false;
		}

		strcpy( inputPQRorPDBFileName, argv[3] );
		strcpy( inputRawSurfaceFileName, argv[4] );
		strcpy( outputMeanRawSurfaceFileName, argv[5] );
		strcpy( outputGaussianRawSurfaceFileName, argv[6] );
		strcpy( curvatureFileName, argv[7] );

		dim1 = atoi(argv[8]); dim2 = atoi(argv[9]); dim3 = atoi(argv[10]);
		blob = atof(argv[11]);
		numberOfGridDivisions = atoi( argv[12] );
		maxFunctionError = atof( argv[13] );
		if( !PDBParser::GroupOfAtoms::intToRadiusType( &radiusType, atoi(argv[14]) ) ) return false;
		level = atoi(argv[15]); 

		return BlurMapsDataManager::getCurvaturesFromIsocontourFile( inputPQRorPDBFileName,
				    				   dim1, dim2, dim3, blob, inputRawSurfaceFileName, 
									   outputMeanRawSurfaceFileName, outputGaussianRawSurfaceFileName,
									   curvatureFileName,
									   numberOfGridDivisions, maxFunctionError, radiusType, level );
	}
	//// new with isovalue
	else if( createIsosurface == 1 )
	{
		if( argc != 16 )
		{
			printUsage();
			return false;
		}
		strcpy( inputPQRorPDBFileName, argv[3] );
		strcpy( outputMeanRawSurfaceFileName, argv[4] );
		strcpy( outputGaussianRawSurfaceFileName, argv[5] );
		strcpy( curvatureFileName, argv[6] );

		dim1 = atoi(argv[7]); dim2 = atoi(argv[8]); dim3 = atoi(argv[9]);
		blob = atof(argv[10]);
		isovalue = atof(argv[11]);
		numberOfGridDivisions = atoi( argv[12] );
		maxFunctionError = atof( argv[13] );
		if( !PDBParser::GroupOfAtoms::intToRadiusType( &radiusType, atoi(argv[14]) ) ) return false;

		level = atoi(argv[15]); 

		SimpleVolumeData* sData = 0;

		///////// create volume from blurring code ////////////
		sData = BlurMapsDataManager::getVolume( inputPQRorPDBFileName, 
			outputVolumeFileName, 
			dim1, dim2, dim3, PDBParser::GroupOfAtoms::ELECTRON_DENSITY,
		    false, blob, PDBParser::GroupOfAtoms::ATOM, NULL, 0, radiusType, level );
		if( !sData ) return false;
		////////////////////////////////////////////////////////


		////// extract an isocontour //////////////////
		Geometry* geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);
		if( geometry == 0 ) return false;
		/////////////////////////////////////////////


		
		/////  get curvatures //////////////
		return BlurMapsDataManager::getCurvatures( inputPQRorPDBFileName, dim1, dim2, dim3, 
			blob, geometry, outputMeanRawSurfaceFileName, 
			outputGaussianRawSurfaceFileName, curvatureFileName,
			numberOfGridDivisions, maxFunctionError, radiusType, level );
		////////////////////////////////////


		return false;
	}
	else
	{
		printUsage();
		return false;
	}
}

//SKVINAY
bool Server::classifyPoints( int argc, char* argv[] )
{
	if( argc != 3 )
	{
		printUsage();
		return false;
	}

	char inputPtsFile[256];
	char outputPtsFile[256];

	strcpy( inputPtsFile, argv[2] );
	strcpy( outputPtsFile, argv[3] );

	bool ret = false;
/*	VorocompDataManager *v = new VorocompDataManager();
	ret = v->classifyPoints( inputPtsFile, outputPtsFile );
	delete v;
*/
	return ret;
}

bool Server::growOut( int argc, char* argv[] )
{
	if( argc != 5 )
	{
		printUsage();
		return false;
	}
	
	char fileNameIn[256];
	char fileNameOut[256];
	double extent;

	strcpy( fileNameIn, argv[2] );
	strcpy( fileNameOut, argv[3] );
	extent = atof(argv[4]);

	DockingManager::growOut( fileNameIn, fileNameOut, extent );
	return false;
}

bool Server::writePDB( int argc, char* argv[] )
{
	if( argc != 5 )
	{
		printUsage();
		return false;
	}

	int outputLevel = atoi(argv[4]);

	///////// create the molecule //////////
	GOALoader* gLoader = new GOALoader();
	PDBParser::GroupOfAtoms* molecule = gLoader->loadFile( argv[2] );
	delete gLoader;
	if( !molecule ) return false;
	////////////////////////////////////////


	///////// write the output ///////////////
	FILE* stream = fopen( argv[3], "w" );
	if( stream == 0 ) 
	{
		delete molecule; molecule= 0;
		return false;
	}
	bool ret = PDBParser::writeGOA2PDB ( stream, molecule, outputLevel, 0);
	fclose( stream );
	//////////////////////////////////////////


	delete molecule; molecule = 0;
	return ret;
}

bool Server::writeGOA( int argc, char* argv[] )
{
	if( argc != 4 )
	{
		printUsage();
		return false;
	}

	///////// create the molecule //////////
	GOALoader* gLoader = new GOALoader();
	PDBParser::GroupOfAtoms* molecule = gLoader->loadFile( argv[2] );
	delete gLoader;
	if( !molecule ) return false;
	////////////////////////////////////////


	///////// write the output ///////////////
	FILE* stream = fopen( argv[3], "w" );
	if( stream == 0 ) 
	{
		delete molecule; molecule= 0;
		return false;
	}
	bool ret = PDBParser::writeGOAtoFile ( stream, molecule);
	fclose( stream );
	//////////////////////////////////////////

	delete molecule; molecule = 0;
	return ret;
}

bool Server::addNewDataSet( int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 3 ) return false;

	char inputFileName[256];

	strcpy( inputFileName, argv[2] );
	return mWindow->addNewDataSet( inputFileName);
}

bool Server::splitView(int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 2 ) return false;

	mWindow->splitViewSlot();
	return true;
}

bool Server::deleteData(int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 3 ) return false;

	int dataSetIndex = atoi(argv[2]);

	mWindow->deleteData( dataSetIndex );
	return true;
}

bool Server::deletePrevData(int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 2 ) return false;

	return mWindow->deletePrevData();
}

bool Server::deleteAllData(int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 2 ) return false;

	return mWindow->deleteAllData();
}

bool Server::setVisible(int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 4 ) return false;
	bool render;
	char renderC[256];

	strcpy( renderC, argv[2] );
#ifdef _WIN32
        if( strcmpi( renderC, "true" ) == 0 )
#else
        if( strcasecmp( renderC, "true" ) == 0 ) 
#endif
		render = true;
	else
		render = false;
	int dataSetIndex = atoi(argv[3]);

	return mWindow->setVisible(render, dataSetIndex);
}

bool Server::setVisiblePrev(int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 3 ) return false;
	bool render;
	char renderC[256];

	strcpy( renderC, argv[2] );
#ifdef _WIN32
        if( strcmpi( renderC, "true" ) == 0 )
#else
        if( strcasecmp( renderC, "true" ) == 0 ) 
#endif
		render = true;
	else
		render = false;

	return mWindow->setVisiblePrev(render);
}

bool Server::saveImage(int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 4 ) return false;
	return mWindow->saveImage(argv[2], argv[3]);
}

bool Server::printPDBInformation( int argc, char* argv[] )
{
	if( argc != 5 ) return false;
	
	BallAndStickDataManager b;
	int printType = atoi(argv[4]);
	return b.printPDBInformation( argv[2], argv[3], printType ); // input , output file names
}

bool Server::setGridVisible(int argc, char* argv[], MoleculeVizMainWindow *mWindow )
{
	if( !mWindow ) return false;
	if( argc != 3 ) return false;

#ifdef _WIN32
        if( strcmpi( argv[2], "true" ) == 0 )
#else
        if( strcasecmp( argv[2], "true" ) == 0 ) 
#endif
		return mWindow->setGridVisible(true);
	else
		return mWindow->setGridVisible(false);
}

bool Server::getMaxDistanceFromPoint( int argc, char* argv[] )
{
	if( argc != 8 ) return false;
	double xOrigin = 0, yOrigin = 0, zOrigin = 0;
	bool appendToFile = true;

	// argv[2], argv[3] = input , output file names
	xOrigin = atof( argv[4] );
	yOrigin = atof( argv[5] );
	zOrigin = atof( argv[6] );

	BallAndStickDataManager b;
#ifdef _WIN32
        if( strcmpi( argv[7], "true" ) == 0 )
#else
        if( strcasecmp( argv[7], "true" ) == 0 ) 
#endif
		return b.getMaxDistanceFromPoint( argv[2], argv[3], xOrigin, yOrigin, zOrigin, true );
	else
		return b.getMaxDistanceFromPoint( argv[2], argv[3], xOrigin, yOrigin, zOrigin, false );
}

bool Server::getHydrophobicityOnSurface( int argc, char* argv[] )
{
	if( argc != 11 ) return false;

	char inputPDBFile[256];
	char inputSurfaceFile[256];
	int dim1,dim2,dim3;
	double blobbiness;
	char outputHydrophobicityValues[256];
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	int level;

	strcpy( inputPDBFile, argv[2] );
	strcpy( inputSurfaceFile, argv[3] );
	dim1 = atoi(argv[4]); dim2 = atoi(argv[5]); dim3 = atoi(argv[6]);
	blobbiness = atof(argv[7]);
	strcpy( outputHydrophobicityValues, argv[8] );
	if( !PDBParser::GroupOfAtoms::intToRadiusType( &radiusType, atoi(argv[9]) ) ) return false;
	level = atoi(argv[10]);
	
	/////// make the hydrophobicity volume //////
	SimpleVolumeData* sData = 0;
	{
		sData = BlurMapsDataManager::getVolume( inputPDBFile, "", 
			dim1, dim2, dim3, PDBParser::GroupOfAtoms::PER_ATOM_HYDROPHOBICITY,
				false, blobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType, level);
	}
	/////////////////////////////////////////////


	/////// read in the isocontour /////////
	SurfaceData* surfaceData = new SurfaceData(0);
	{
		surfaceData->read( inputSurfaceFile );
	}
	////////////////////////////////////////


	//////// error checking ///////////////
	{
		if( !sData || !surfaceData )
		{
			delete sData; sData = 0;
			delete surfaceData; surfaceData = 0;
			return false;
		}
		if( !surfaceData->getGeometry() )
		{
			delete sData; sData = 0;
			delete surfaceData; surfaceData = 0;
			return false;
		}
		if( surfaceData->getGeometry()->m_NumTriVerts < 1 )
		{
			delete sData; sData = 0;
			delete surfaceData; surfaceData = 0;
			return false;
		}
	}
	////////////////////////////////////////


	////// get the hydrophobicity values for each vertex ////////
	double* funcVals = new double[surfaceData->getGeometry()->m_NumTriVerts];
	{
		if( !SimpleVolumeDataIsocontourer::getFunctionValues( sData, surfaceData->getGeometry(), funcVals ) )
		{
			delete funcVals;
			delete sData; sData = 0;
			delete surfaceData; surfaceData = 0;
			return false;
		}
	}
	/////////////////////////////////////////////////////////////


	////// output the values /////////////
	{
		FILE* fp = fopen( outputHydrophobicityValues, "w" );
		if( !fp )
		{
			delete funcVals;
			delete sData; sData = 0;
			delete surfaceData; surfaceData = 0;
			return false;
		}
		fprintf( fp, "%d\n", surfaceData->getGeometry()->m_NumTriVerts );
		int i;
		for( i=0; i<surfaceData->getGeometry()->m_NumTriVerts; i++ )
		{
			fprintf( fp, "%lf\n", funcVals[i] );
		}
		fclose( fp );
	}
	//////////////////////////////////////

	///// assign hydrophobicity values to the surface ///////////
	{
		if( surfaceData->getGeometry()->m_NumTriVerts < 1 ) return false;

		float r, g, b;
		int i;

		delete []surfaceData->getGeometry()->m_TriVertColors;
		surfaceData->getGeometry()->m_TriVertColors = new float[surfaceData->getGeometry()->m_NumTriVerts*3];

		for( i=0; i<surfaceData->getGeometry()->m_NumTriVerts; i++ )
		{	
			r = g = b = 1.0;
			double v = funcVals[i];
			if( v < 0 ) 
			{
				g = 1.0 + v*5;
				b = 1.0 + v*5;
				if( v < -0.2 )
				{
					g = 0.0; b = 0.0; 
				}
			}
			if( v > 0 ) 
			{
				r = 1.0 - v*5;
				b = 1.0 - v*5;
				if( v > 0.2 )
				{
					r = 0.0; b = 0.0; 
				}
			}
			surfaceData->getGeometry()->m_TriVertColors[3*i+0] = r;
			surfaceData->getGeometry()->m_TriVertColors[3*i+1] = g;
			surfaceData->getGeometry()->m_TriVertColors[3*i+2] = b;
		}

		GeometryLoader* geometryLoader = new GeometryLoader();
		bool savedMeanCurvatureFile = geometryLoader->saveFile( inputSurfaceFile, "Rawnc files (*.rawnc)", surfaceData->getGeometry() );
	}

	delete surfaceData;
	delete funcVals;
	delete sData;

	return true;
}

bool Server::getSurfaceFromPDB( int argc, char* argv[] )
{
	if( argc != 11 )
	{
		printUsage();
		return false;
	}

	char pdbFileName[256];
	char surfaceFileName[256];
	double isovalue = 1;
	int dim1, dim2, dim3;
	double blobbiness;
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	int level;

	strcpy( pdbFileName, argv[2] );
	strcpy( surfaceFileName, argv[3] );
	isovalue = atof( argv[4] );
	dim1 = atoi(argv[5]); dim2 = atoi(argv[6]); dim3 = atoi(argv[7]);
	blobbiness = atof(argv[8]);
	if( !PDBParser::GroupOfAtoms::intToRadiusType( &radiusType, atoi(argv[9]) ) ) return false;
	level = atoi(argv[10]);

	///////// create the volume ////////////////////
	SimpleVolumeData* sData = 0;
	{
		sData = BlurMapsDataManager::getVolume( pdbFileName, "", 
			    dim1, dim2, dim3, PDBParser::GroupOfAtoms::ELECTRON_DENSITY,
				false, blobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType, level);
		if( !sData )
			return false;
	}
	/////////////////////////////////////////////////


	////// create the isosurface /////////////////////
	Geometry* geometry = 0;
	{
		geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);
		if( geometry == 0 ) 
		{
			delete sData;
			return false;
		}
	}
	///////////////////////////////////////////////////


	///////// save the isocontour /////////////////////
	{
		GeometryLoader* geometryLoader = new GeometryLoader();
		if( !geometryLoader->saveFile( surfaceFileName, "Rawnc files (*.rawnc)", geometry) )
		{
			delete geometry;
			delete sData;
			return false;
		}
		delete geometryLoader;
	}
	///////////////////////////////////////////////////


	delete geometry;
	delete sData;
	return true;
}

bool Server::createSkinRegion( int argc, char* argv[] )
{
	if( argc != 9 ) 
	{
		printUsage();
		return false;
	}

	char inputFileName[256];
	char outputFileName[256];
	int dim1, dim2, dim3;
	double probeRadius;
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;
	int depth;

	strcpy( inputFileName, argv[2] );
	strcpy( outputFileName, argv[3] );
	dim1 = atoi(argv[4]); dim2 = atoi(argv[5]); dim3 = atoi(argv[6]);
	probeRadius = atof(argv[7]);
	if( !PDBParser::GroupOfAtoms::intToRadiusType( &radiusType, atoi(argv[8]) ) ) return false;

	depth = 2;
	SimpleVolumeData* sData;

	sData = BlurMapsDataManager::getSkinRegionVolume( inputFileName, outputFileName, 
		    dim1, dim2, dim3, probeRadius, radiusType, depth);

	bool ret;
	if( sData ) 
		ret = true;
	else
		ret = false;

	delete sData;

	return ret;
}

bool Server::correlate( int argc, char* argv[] )
{
	if( argc != 10 ) 
	{
		printUsage();
		return false;
	}

	char hydroFileName[256];
	char curvFileName[256];

	bool discretizeX;
	bool discretizeY;

	double discreteMaxVal;
	double discreteMinVal;

	int rangeToCorrelate; // 0 -> use all, < 0 -> use neg, > 0 -> use pos

	char outFileName[256];

	strcpy( hydroFileName, argv[2] );
	strcpy( curvFileName, argv[3] );

	if( strcmp( argv[4], "true" ) == 0 ) 
		discretizeX = true;
	else
		discretizeX = false;

	if( strcmp( argv[5], "true" ) == 0 ) 
		discretizeY = true;
	else
		discretizeY = false;

	discreteMaxVal = atof(argv[6]);
	discreteMinVal = atof(argv[7]);

	rangeToCorrelate = atoi(argv[8]);

	strcpy( outFileName, argv[9] );

	// read allocate arrays
	double* h = 0;
	double* m = 0;
	double* g = 0;
	int  n = 0;

	{
		FILE* fp1 = fopen( hydroFileName, "r" );
		if( !fp1 ) return false;
		if( !fscanf( fp1, "%d\n", &n ) == 1) return false;
		if( n < 1 ) return false;

		FILE* fp2 = fopen( curvFileName, "r" );
		if( !fp2 ) return false;
		int temp;
		if( !fscanf( fp2, "%d\n", &temp ) == 1) return false;
		if( temp != n ) return false;

		h = new double[n]; m = new double[n]; g = new double[n];
		int i;
		for( i=0; i<n; i++ )
		{
			if( fscanf( fp1, "%lf\n", &(h[i]) ) != 1 ) return false;
			if( fscanf( fp2, "%lf %lf\n", &(m[i]), &(g[i]) ) != 2 ) return false;
		}
		fclose( fp1 );
		fclose( fp2 );
	}

	// discretize if required
	if( discretizeX )
	{
		if( !CCVOpenGLMath::LinearAlgebra::discretize( h, n, discreteMaxVal, discreteMinVal ) ) return false;
	}
	if( discretizeY )
	{
		if( !CCVOpenGLMath::LinearAlgebra::discretize( m, n, discreteMaxVal, discreteMinVal ) ) return false;
		if( !CCVOpenGLMath::LinearAlgebra::discretize( g, n, discreteMaxVal, discreteMinVal ) ) return false;
	}

	// correllate
	double coef_hydro_mean = 0;
	double coef_hydro_gaus = 0;

	if( !CCVOpenGLMath::LinearAlgebra::selectivelyCorrelate( m, rangeToCorrelate, h, n, &coef_hydro_mean ) ) return false;
	if( !CCVOpenGLMath::LinearAlgebra::selectivelyCorrelate( g, rangeToCorrelate, h, n, &coef_hydro_gaus ) ) return false;

	// print out result
	{
		FILE* fp = fopen( outFileName,"a");
		fprintf(fp, "---------------------------------------\n");
		fprintf(fp, "Correllating:\n");
		fprintf(fp, "\t%s %s\n", hydroFileName, curvFileName );
		if( discretizeX )
			fprintf(fp, "\tDiscretized hydro with values %lf %lf\n", discreteMaxVal, discreteMinVal );
		if( discretizeY )
			fprintf(fp, "\tDiscretized curvs with values %lf %lf\n", discreteMaxVal, discreteMinVal );
		if( rangeToCorrelate == 0 )
			fprintf(fp, "\tCorrellated the entire array\n");
		if( rangeToCorrelate < 0 )
			fprintf(fp, "\tCorrellated the negative hydro only\n");
		if( rangeToCorrelate > 0 )
			fprintf(fp, "\tCorrelated the positive hydro only\n");
		fprintf(fp, "\tDot product <hydro.mean>,<hydro.gaus> = %lf %lf\n\n", coef_hydro_mean, coef_hydro_gaus );
		fclose( fp );
	}
	return true;
}

bool Server::getArea( int argc, char* argv[] )
{
	if( argc != 3 ) 
	{
		printUsage();
		return false;
	}


	/////// read in the isocontour /////////
	SurfaceData* surfaceData = new SurfaceData(0);
	{
		if( !surfaceData->read( argv[2] ) )
		{
			delete surfaceData; surfaceData = 0;
			printf("Could not read file %s\n", argv[2] );
			return false;
		}
	}
	////////////////////////////////////////


	/////// calculate surface area and print it out ////////
	printf("Area of surface %s = %lf\n", argv[2], surfaceData->getArea() );
	delete surfaceData;
	////////////////////////////////////////////////////////

	return true;
}

bool Server::getPatchAreas( int argc, char* argv[] )
{
	if( argc != 6 ) 
	{
		printUsage();
		return false;
	}

	double isovalue;

	isovalue = atof( argv[4] );

	/////// read in the isocontour /////////
	SurfaceData* surfaceData = new SurfaceData(0);
	{
		if( !surfaceData->read( argv[2] ) )
		{
			delete surfaceData; surfaceData = 0;
			printf("Could not read file %s\n", argv[2] );
			return false;
		}
		if( !surfaceData->getGeometry() )
		{
			delete surfaceData; surfaceData = 0;
			printf("Could not read file %s\n", argv[2] );
			return false;
		}
	}
	////////////////////////////////////////


	/////// read the function values /////////
	double* functionAtVertices = 0;
	{
		FILE* fp = fopen( argv[3], "r" );
		if( !fp )
		{
			delete surfaceData; surfaceData = 0;
			printf("Could not read file %s\n", argv[3] );
			return false;
		}
		int numVerts = 0;
		fscanf( fp, "%d\n", &numVerts );
		if( (numVerts < 1) || (numVerts!=surfaceData->getGeometry()->m_NumTriVerts) )
		{
			delete surfaceData; surfaceData = 0;
			printf("Number of vertices in file %s was wrong\n", argv[3] );
			return false;
		}

		functionAtVertices = new double[numVerts];
		int i;
		for( i=0; i<numVerts; i++ )
		{
			if( fscanf( fp, "%lf\n", &(functionAtVertices[i]) ) != 1 )
			{
				delete surfaceData; surfaceData = 0;
				delete []functionAtVertices; functionAtVertices = 0;
				printf("Line %d in file %s was wrong\n", i+1, argv[3] );
				return false;
			}
		}

		fclose( fp );
	}
	//////////////////////////////////////////


	/////// calculate surface areas and print it out ////////
	double areaBelow = 0;
	double areaAbove = 0;
	if( !surfaceData->getAreas(isovalue, functionAtVertices, &areaBelow, &areaAbove ) )
	{
		printf("Could not compute area of surface %s\n", argv[2] );
		delete surfaceData; surfaceData = 0;
		delete []functionAtVertices; functionAtVertices = 0;
		return false;
	}
	{
		FILE* fp = fopen( argv[5], "a" );
		if( !fp )
		{
			printf("Could not open file%s\n", argv[5] );
			delete surfaceData; surfaceData = 0;
			delete []functionAtVertices; functionAtVertices = 0;
			return false;
		}

		fprintf( fp, "Surface: %s,   Isovalue: %lf,  Area below, above = %lf, %lf\n", argv[2], isovalue, areaBelow, areaAbove );
		fclose( fp );
	}
	/////////////////////////////////////////////////////////


	delete []functionAtVertices; functionAtVertices = 0;
	delete surfaceData; surfaceData = 0;
	return true;
}

bool Server::populateSAS( int argc, char* argv[] )
{
	if( argc != 9 ) 
	{
		printUsage();
		return false;
	}

	char inputFileName[256];
	char outputFileName[256];
	int dim1, dim2, dim3;
	double probeRadius;
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType;

	strcpy( inputFileName, argv[2] );
	strcpy( outputFileName, argv[3] );
	dim1 = atoi(argv[4]); dim2 = atoi(argv[5]); dim3 = atoi(argv[6]);
	probeRadius = atof(argv[7]);
	if( !PDBParser::GroupOfAtoms::intToRadiusType( &radiusType, atoi(argv[8]) ) ) return false;

	return BlurMapsDataManager::populateSAS( inputFileName, outputFileName, dim1, dim2, dim3, probeRadius, radiusType);
}

bool Server::addVolumes( int argc, char* argv[] )
{
	if( argc != 7 ) 
	{
		printUsage();
		return false;
	}

	char inputFileName1[256];
	char inputFileName2[256];
	char outputFileName[256];
	double scale;
	double sum;

	strcpy( inputFileName1, argv[2] );
	strcpy( inputFileName2, argv[3] );
	strcpy( outputFileName, argv[4] );
	scale = atof(argv[5]);
	sum = atof(argv[6]);

	VolumeLoader* volumeLoader1 = new VolumeLoader();
	SimpleVolumeData* sData1 = volumeLoader1->loadFile( inputFileName1 );
	delete volumeLoader1;
	if( !sData1 ) return false;

	VolumeLoader* volumeLoader2 = new VolumeLoader();
	SimpleVolumeData* sData2 = volumeLoader2->loadFile( inputFileName2 );
	delete volumeLoader2;
	if( !sData2 ) 
	{
		delete sData1;
		return false;
	}

	if( !sData1->addVolume( sData2, scale, sum ) )
	{
		delete sData1;
		delete sData2;
		return false;
	}

	VolumeLoader* volumeLoader = new VolumeLoader();
	bool ret = volumeLoader->saveFile( outputFileName, sData1 );
	delete volumeLoader;
	delete sData1;
	delete sData2;

	return ret;
}

bool getExtension( char* filename, char* extension )
{
	int l = strlen( filename );
	int i;
	for( i=l-1; i>=0; i-- )
	{
		if( filename[i] == '.' ) break;
	}
	if( i == -1 ) return false;
	if( i == l-1 ) return false;

	{
		int j;
		for( j=i+1; j<l; j++ )
			extension[j-(i+1)] = filename[j];
		extension[j-(i+1)] =  '\0';
	}
	return true;
}

bool Server::convert( int argc, char* argv[] )
{
	if( argc < 4 ) 
	{
		printUsage();
		return false;
	}

	char inputFileName[256];
	char outputFileName[256];
	char inputFileNameExtension[256];
	char outputFileNameExtension[256];

	strcpy( inputFileName, argv[2] );
	strcpy( outputFileName, argv[3] );
	if( !getExtension( inputFileName, inputFileNameExtension ) ) return false;
	if( !getExtension( outputFileName, outputFileNameExtension ) ) return false;


	{
		VolumeLoader* vLoader = new VolumeLoader();
		if( vLoader->isValidExtension( inputFileNameExtension ) && vLoader->isValidExtension( outputFileNameExtension ) )
		{
			bool ret = vLoader->saveFile( outputFileName, vLoader->loadFile( inputFileName ) );
			delete vLoader;
			return ret;
		}
		delete vLoader;
	}
	{
		GeometryLoader* gLoader = new GeometryLoader();
		if( gLoader->isValidExtension( inputFileNameExtension ) && gLoader->isValidExtension( outputFileNameExtension ) )
		{
			if( argc == 4 )
			{
				bool ret = gLoader->saveFile( outputFileName, gLoader->loadFile( inputFileName ) );
				delete gLoader;
				return ret;
			}
			else
			{
				Geometry* geometry = gLoader->loadFile( inputFileName );
				SurfaceData* surfaceData = new SurfaceData(0);
				surfaceData->readTransformations( argv[4] );
				surfaceData->setGeometry(geometry);
				Geometry* dup_geometry = surfaceData->createDuplicates();
				delete surfaceData;
				bool ret = gLoader->saveFile( outputFileName, dup_geometry );
				delete gLoader;
				return ret;
			}
		}
		delete gLoader;
	}
	{
		GOALoader* gLoader = new GOALoader();
		if( gLoader->isValidExtension( inputFileNameExtension ) && gLoader->isValidExtension( outputFileNameExtension ) )
		{
			bool ret = gLoader->saveFile( outputFileName, gLoader->loadFile( inputFileName ), PDBParser::ATOM_TYPE, 0 );
			delete gLoader;
			return ret;
		}
		delete gLoader;
	}

	return false;
}

bool isClose( double x1, double y1, double z1, double r, double x2, double y2, double z2, double* d )
{
	double d1 = (x1-x2)*(x1-x2) + 
				(y1-y2)*(y1-y2) + 
				(z1-z2)*(z1-z2);
	*d = sqrt(d1);
	if( d1 < r*r  )
	{
		//printf( "dist = %lf r = %lf\n", sqrt(d1), r );
		return true;
	}
	return false;
}

/*
	Creates:
	1. A file as input for the docking code.
	2. A pdb with a corresponding ".bdy" file.
*/
bool Server::getSurfaceAtoms( int argc, char* argv[] )
{
	if( argc != 4 ) 
	{
		printUsage();
		return false;
	}

	char inputFileName[256];
	char outputFileName[256];

	strcpy( inputFileName, argv[2] );
	strcpy( outputFileName, argv[3] );

	GOALoader* gLoader = new GOALoader();
	PDBParser::GroupOfAtoms* molecule = gLoader->loadFile( inputFileName );
	delete gLoader;

	if( !molecule ) return false;

	double probeRadius = 1.4;
	PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType = PDBParser::GroupOfAtoms::VDW_RADIUS;
	bool* boundaryAtom = 0;
	int depth = 1;

	int dim1 = 128;
	int dim2 = 128;
	int dim3 = 128;
	SimpleVolumeData* sData = new SimpleVolumeData(dim1, dim2, dim3);
	SkinRegion2* skin = new SkinRegion2();
	if( !skin->getSkinRegion( molecule, dim1, dim2, dim3, probeRadius, radiusType, sData, depth ) )
	{
		delete sData;
		delete skin;
		delete molecule;
		return false;
	}
	delete skin;

	Geometry* geometry = 0;
	{
		float isovalue = 1.2f;
		geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, isovalue);
		if( !geometry ) 
		{
			delete molecule;
			delete sData;
			return false;
		}
	}
	delete sData;

	// the code below breaks up isocontour into separate components
	vector<Geometry *> components;

	geometry->separateComponents( &components );
	delete geometry;
	if( components.size() < 1 )
	{
		int i;
		for( i=0; i<components.size(); i++ )
			delete components[i];
		components.clear();
		delete molecule;
		return false;
	}

	Geometry* largestComponent = 0;
	{
		int curSize = components[0]->m_NumTris;
		largestComponent = components[0];
		int i;
		for( i=1; i<components.size(); i++ )
		{
			if( components[i]->m_NumTris > curSize )
			{
				curSize = components[i]->m_NumTris;
				largestComponent = components[i];
			}
		}
		for( i=1; i<components.size(); i++ )
		{
			if( components[i] != largestComponent )
				delete components[i];
		}
		components.clear();
	}

	// find surface atoms of molecule as those close to the largest component
	vector<PDBParser::Atom *> atomList;
	{
		PDBParser::CollectionData* collectionData = 0;
		if( molecule->type == PDBParser::COLLECTION_TYPE ) 
			collectionData = molecule->m_CollectionData;
		BlurMapsDataManager::flattenGOA(molecule, atomList, collectionData, 0, 0, 0, radiusType, PDBParser::ATOM_TYPE, false  );
	}
	boundaryAtom = new bool[atomList.size()];

	{
		int i;
		for( i=0; i<atomList.size(); i++ )
		{
			boundaryAtom[i] = false;

			int v;
			for( v=0; v<largestComponent->m_NumTriVerts; v++ )
			{
				double dist = 10;
				if( isClose( atomList[i]->m_Position[0],
							 atomList[i]->m_Position[1],
							 atomList[i]->m_Position[2],
							 atomList[i]->getRadius() + 0.2,
							 largestComponent->m_TriVerts[v*3+0],
							 largestComponent->m_TriVerts[v*3+1],
							 largestComponent->m_TriVerts[v*3+2],
							 &dist
							 ) )
				{
					//printf("Atom %d is %lf away from surface\n", i, dist );
					boundaryAtom[i] = true;
					break;
				}
			}
		}
	}
	delete molecule;
	delete largestComponent;

	{
		int numBdyAtoms = 0;
		
		// the boundary atoms file takes in surface area exposed for the given atom.
		char bdyFile[256];
		strcpy( bdyFile, inputFileName );
		strcat( bdyFile, ".bdy" );
		FILE* fp = fopen( bdyFile, "w" );

		char interiorfname[256];
		strcpy( interiorfname, outputFileName );
		strcat( interiorfname, "_interior.pdb");
		FILE* fpOut = fopen( interiorfname, "w" );

		char skinfname[256];
		strcpy( skinfname, outputFileName );
		strcat( skinfname, "_skin.pdb");
		FILE* fpSkin = fopen( skinfname,"w");

		int i;
		for( i=0; i<atomList.size(); i++ )
		{
			if( boundaryAtom[i] )
			{
				fprintf( fp, "100000.0\n");
				fprintf(fpSkin, "ATOM  11111  CA  GLU A1111    %8.3f%8.3f%8.3f%6.2f                \n", atomList[i]->m_Position[0], atomList[i]->m_Position[1], atomList[i]->m_Position[2], atomList[i]->getRadius() );
				numBdyAtoms++;
			}
			else
			{
				fprintf( fp, "0.0\n");
				fprintf(fpOut, "ATOM  11111  N   GLU A1111    %8.3f%8.3f%8.3f%6.2f                \n", atomList[i]->m_Position[0], atomList[i]->m_Position[1], atomList[i]->m_Position[2], atomList[i]->getRadius() );	
			}
		}
		fclose( fp );
		fclose( fpOut );
		fclose( fpSkin );

		{
			// write docking input file
			FILE* fp = fopen( outputFileName, "w" );
			if( !fp  )return false;
			
			fprintf(fp, "%d\n", atomList.size() );
			int i;
			for( i=0; i<atomList.size(); i++ )
			{
				if( boundaryAtom[i] ) continue;

				PDBParser::Atom *at = atomList[i];

				fprintf( fp, "I     %12.5lf %12.5lf %12.5lf %12.5f\n", at->m_Position[0], at->m_Position[1], at->m_Position[2], at->getCharge() );
			}

			for( i=0; i<atomList.size(); i++ )
			{
				if( boundaryAtom[i] ) 
				{
					PDBParser::Atom *at = atomList[i];

					fprintf( fp, "E     %12.5lf %12.5lf %12.5lf %12.5f\n", at->m_Position[0], at->m_Position[1], at->m_Position[2], at->getCharge() );
				}
			}
			fclose( fp );
		}
	}

	{
		int n = atomList.size();
		int i;
		for( i=0; i<n; i++ )
			delete atomList[i];
		atomList.clear();
	}
	if( boundaryAtom )
	{
		delete [] boundaryAtom; boundaryAtom = 0;
	}
	return true;
}

/*void Server::getContourSpectrum(unsigned char* uchar_data, int type, int* dim, int array_size, float* isoval , float* area, float* min_vol, float* max_vol, float* gradient)
{
	int i;
	ConDataset* the_data;
	the_data = newDatasetReg(type, CONTOUR_REG_3D, 1, 1, dim,uchar_data);
	
	Signature       *sig;
	sig=getSignatureFunctions(the_data, 0,0);
	
	for (i=0;i<array_size;i++) {
		isoval[i]=sig[0].fx[i];
		area[i]=sig[0].fy[i];
		min_vol[i]=sig[1].fy[i];
		max_vol[i]=sig[2].fy[i];
		gradient[i]=sig[3].fy[i];
	}
	
	delete the_data;
	delete sig;
}*/

bool Server::getContourStats(int argc, char* argv[] )
{
	if( argc != 4 && argc != 5) 
	{
		printUsage();
		return false;
	}

	char inputFileName[256];
	char outputFileName[256];

	strcpy( inputFileName, argv[2] );
	strcpy( outputFileName, argv[3] );

	SimpleVolumeData* sData = 0;
	double isovalue = 1;

	if( argc == 4 )
	{
		bool useBlurring = true;
		///////// create the volume ////////////////////
		int dim1 = 128, dim2 = 128, dim3 = 128;
		if( useBlurring )
		{
			float blobbiness = -2.3f;
			int level = PDBParser::ATOM_TYPE ;
			PDBParser::GroupOfAtoms::RADIUS_TYPE radiusType = PDBParser::GroupOfAtoms::VDW_RADIUS;
			sData = BlurMapsDataManager::getVolume( inputFileName, "", 
					dim1, dim2, dim3, PDBParser::GroupOfAtoms::ELECTRON_DENSITY,
					false, blobbiness, PDBParser::GroupOfAtoms::ATOM, 0, 0, radiusType, level);
			if( !sData )
				return false;
		}
		else
		{
			double probeRadius = 1.4;
			isovalue = probeRadius;
			int depth = 3;

			GOALoader* gLoader = new GOALoader();
			PDBParser::GroupOfAtoms* molecule = gLoader->loadFile( inputFileName );
			delete gLoader;
			if( !molecule ) return false;

			sData = new SimpleVolumeData(dim1, dim2, dim3);
			SkinRegion2* skin = new SkinRegion2();
			if( !skin->getSkinRegion( molecule, dim1, dim2, dim3, probeRadius, PDBParser::GroupOfAtoms::VDW_RADIUS, sData, depth ) )
			{
				delete sData; sData = 0;
				delete skin; skin = 0;
				return false;
			}
			delete skin;	
		}
		/////////////////////////////////////////////////
	}
	else
	{
		VolumeLoader* vLoader = new VolumeLoader();
		sData = vLoader->loadFile( inputFileName );
		delete vLoader;

		if( !sData )
			return false;

		isovalue = atof( argv[4]);
	}

	/////////// get the contour spectrum //////////
	if( sData->getNumberOfVariables() < 0 ) return false;

	int type;
	void* data;
	
	double minVal = 0, maxVal = 1;
	if( sData->getNumberOfVariables() < 4 ) 
	{
		// rawiv
		data = sData->getData(0);
		type = sData->getType(0);
		if( !sData ->getMinMax( 0, &minVal, &maxVal ) ) return false;
	}
	else
	{
		// rawv
		data = sData->getData(3);
		type = sData->getType(3);
		if( !sData->getMinMax( 3, &minVal, &maxVal ) ) return false;
	}
	
	int dim[3];

	dim[0] = sData->getWidth();
	dim[1] = sData->getHeight();
	dim[2] = sData->getDepth();

	float span[3];
	span[0] = sData->getSpanX();
	span[1] = sData->getSpanY();
	span[2] = sData->getSpanZ();

	float orig[3];
	orig[0] = sData->getMinX();
	orig[1] = sData->getMinY();
	orig[2] = sData->getMinZ();

	float* vol = 0;

	{
		int i;
		vol = new float[dim[0]*dim[1]*dim[2]];

		switch (type) {
		case SimpleVolumeData::UCHAR:
			{
				unsigned char* ucData = (unsigned char*)data;
				for( i=0; i<dim[0]*dim[1]*dim[2]; i++ )
					vol[i] = ucData[i];
			}
			break;
		case SimpleVolumeData::USHORT:
			{
				unsigned short* usData = (unsigned short*)data;
				for( i=0; i<dim[0]*dim[1]*dim[2]; i++ )
					vol[i] = usData[i];
			}
			break;
		case SimpleVolumeData::ULONG:
			{
				unsigned long* ulData = (unsigned long*)data;
				for( i=0; i<dim[0]*dim[1]*dim[2]; i++ )
					vol[i] = ulData[i];
			}
			break;
		case SimpleVolumeData::FLOAT:
			{
				float* fData = (float*)data;
				for( i=0; i<dim[0]*dim[1]*dim[2]; i++ )
					vol[i] = fData[i];
			}
			break;
		case SimpleVolumeData::DOUBLE:
			{
				double* dData = (double*)data;
				for( i=0; i<dim[0]*dim[1]*dim[2]; i++ )
					vol[i] = dData[i];
			}
			break;
			
		default:
			return false;
		}
	}
	delete sData;
	/////////////////////////////////////////////////////////////////

	double area = 0;
	double volume = 0;
	AreaVolume* areaVolume = new AreaVolume();
	areaVolume->getVolume(vol, isovalue, dim, orig, span, &area, &volume);
	delete areaVolume;


	FILE* fp = fopen(outputFileName, "a");
	fprintf(fp, "area, volume = %015lf %015lf\n", area, volume );

	fclose( fp );

	return true;
}

bool Server::downloadPDB(int argc, char* argv[], MoleculeVizMainWindow *mWindow  )
{
	if( argc != 4 ) 
	{
		printUsage();
		return false;
	}

	char pdbID[256];
	char outputFileName[256];

	strcpy( pdbID, argv[2] );
	strcpy( outputFileName, argv[3] );

	DownloadPDB* downloadPDB = new DownloadPDB();

	if( !downloadPDB->blockedDownload( pdbID, outputFileName ) )
	{
		delete downloadPDB;
		return false;
	}

	return true;
}

bool Server::getSignDistanceFunction(int argc, char* argv[] )
{
	if( argc != 6) 
	{
		printUsage();
		return false;
	}

	char inputFileName[256];
	char outputVolumeFileName[256];
	int size = 32;
	int flipNormals = 0;

	strcpy( inputFileName, argv[2] );
	strcpy( outputVolumeFileName, argv[3] );
	size = atoi(argv[4]); 
	flipNormals = atoi(argv[5]); 
	
	char inputFileNameExtension[256];

	if( !getExtension( inputFileName, inputFileNameExtension ) ) return false;

	Geometry* geometry = 0;
	if( !geometry )
	{
		VolumeLoader* vLoader = new VolumeLoader();
		if( vLoader->isValidExtension( inputFileNameExtension ))
		{
			// load volume, isocontour it.
			SimpleVolumeData* sData = vLoader->loadFile( inputFileName );
			if( !sData )return false;
			geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, 1.0);
		}
		delete vLoader;
	}
	if( !geometry )
	{
		GeometryLoader* gLoader = new GeometryLoader();
		if( gLoader->isValidExtension( inputFileNameExtension ) )
		{
			GeometryLoader* gLoader = new GeometryLoader();
			geometry = gLoader->loadFile( inputFileName );
			delete gLoader;
		}
		delete gLoader;
	}
	if( !geometry )
	{
		GOALoader* gLoader = new GOALoader();
		if( gLoader->isValidExtension( inputFileNameExtension ) )
		{
			// load goa. Blur it. Isocontour it.
			SimpleVolumeData* sData = 0;
			sData = BlurMapsDataManager::getVolume( inputFileName, "", 
					64, 64, 64, PDBParser::GroupOfAtoms::ELECTRON_DENSITY,
					false, -2.3, PDBParser::GroupOfAtoms::ATOM, 0, 0, 
					PDBParser::GroupOfAtoms::VDW_RADIUS, 0);
			if( !sData ) return false;
			geometry = SimpleVolumeDataIsocontourer::getIsocontour(sData, 1.0);
			delete sData;
		}
		delete gLoader;
	}

	if( !geometry ) return false;

	float mins[3];
	float maxs[3];
	geometry->CalculateExtents();
	mins[0] = geometry->m_Min[0]; mins[1] = geometry->m_Min[1]; mins[2] = geometry->m_Min[2];
	maxs[0] = geometry->m_Max[0]; maxs[1] = geometry->m_Max[1]; maxs[2] = geometry->m_Max[2];


	SDFLibrary::setParameters(size, flipNormals, mins, maxs);

	//Finally, call the SDF code.
	float* values = SDFLibrary::computeSDF(geometry->m_NumTriVerts, geometry->m_TriVerts, geometry->m_NumTris, (int*)(geometry->m_Tris));
	delete geometry;
	if( !values ) return false;

	float* choppedValues = new float[size*size*size];
	{
		int i, j, k;
		int c=0;
		for( i=0; i<=size; i++ )
		{
			for( j=0; j<=size; j++ )
			{
				for( k=0; k<=size; k++ )
				{
					if( i!=size && j!=size && k!=size )
						choppedValues[c++] = values[i*(size+1)*(size+1) + j*(size+1) + k];
				}
			}
		}
	}
	delete []values;

	//Lastly, the header info.
	SDFLibrary::RAWIV_header* volInfo = SDFLibrary::getVolumeInfo();
	if( !volInfo ) return false;

	SimpleVolumeData* sData = 0;

	{
		unsigned int dims[3]; dims[0] = size; dims[1] = size; dims[2] = size;
		sData = new SimpleVolumeData(dims); // Johns oddities
		sData->setDimensions( dims );
	}
	sData->setNumberOfVariables(1);
	sData->setData(0, choppedValues);
	sData->setType(0, SimpleVolumeData::FLOAT);
	sData->setName(0, "Sign distance function from TexMol");
	sData->setMinExtent(volInfo->minext);
	sData->setMaxExtent(volInfo->maxext);

	VolumeLoader* vLoader = new VolumeLoader();
	bool ret = vLoader->saveFile( outputVolumeFileName, sData );
	delete vLoader;

	delete sData;

	return ret;
}

bool Server::getMolecularCharacteristics(int argc, char* argv[] )
{
	if( argc < 6) 
	{
		printUsage();
		return false;
	}

	char inputGOAFileName[256];
	char inputElecVolumeFileName[256];
	char outputFileName[256];
	int numberOfSpheres = 0;
	float* activeSite = 0;
	double* radii = 0;


	strcpy( inputGOAFileName, argv[2] );
	strcpy( inputElecVolumeFileName, argv[3] );
	strcpy( outputFileName, argv[4] );
	numberOfSpheres = atoi(argv[5]); 

	if( numberOfSpheres > 0 )
	{
		if( numberOfSpheres + 6 + 3 != argc ) return false;

		activeSite = new float[3];
		radii = new double[numberOfSpheres];
	
		activeSite[0] = atof( argv[6] );
		activeSite[1] = atof( argv[7] );
		activeSite[2] = atof( argv[8] );

		int i;
		for( i=0; i<numberOfSpheres; i++ )
			radii[i] = atof( argv[9+i] );
	}

	// note, molecule OR volume can be NULL !
	GOALoader* gLoader = new GOALoader();
	PDBParser::GroupOfAtoms* molecule = gLoader->loadFile( inputGOAFileName );
	delete gLoader;

	VolumeLoader* volumeLoader = new VolumeLoader();
	SimpleVolumeData* electrostatics = volumeLoader->loadFile( inputElecVolumeFileName );
	delete volumeLoader;

	if( !molecule && !electrostatics )
	{
		if( activeSite ) { delete []activeSite; activeSite = 0; }
		if( radii ) { delete []radii; radii = 0; }
		return false;
	}

	MolecularCharacteristics* molecularCharacteristics = new MolecularCharacteristics(	
								molecule, 
								electrostatics,
								numberOfSpheres,
								activeSite,
								radii);

	if( !molecularCharacteristics->appendHeader( outputFileName ) )
	{
		if( activeSite ) { delete []activeSite; activeSite = 0; }
		if( radii ) { delete []radii; radii = 0; }
		return false;
	}

	if( !molecularCharacteristics->appendCharacteristics( outputFileName ) )
	{
		if( activeSite ) { delete []activeSite; activeSite = 0; }
		if( radii ) { delete []radii; radii = 0; }
		return false;
	}

	if( activeSite ) { delete []activeSite; activeSite = 0; }
	if( radii ) { delete []radii; radii = 0; }
	return true;
}

bool Server::expandMolecule( int argc, char* argv[] )
{
	// read a GOA file. Apply the transformations and save into a new GOA file
	if( argc < 4) 
	{
		printUsage();
		return false;
	}

	char inputGOAFileName[256];
	char outputGOAFileName[256];


	strcpy( inputGOAFileName, argv[2] );
	strcpy( outputGOAFileName, argv[3] );

	PDBParser::GroupOfAtoms* molecule = 0;
	/////// read file
	{
		GOALoader* gLoader = new GOALoader();
		molecule = gLoader->loadFile( inputGOAFileName );
		delete gLoader;
		if( !molecule ) return false;
	}

	vector<PDBParser::Atom *> atomList;

	// "flatten" the GOA
	PDBParser::CollectionData* collectionData = 0;
	if( molecule->type == PDBParser::COLLECTION_TYPE ) 
		collectionData = molecule->m_CollectionData;
	BlurMapsDataManager::flattenGOA(molecule, atomList, collectionData, 0, 0, 0, 
		PDBParser::GroupOfAtoms::VDW_RADIUS, PDBParser::GroupOfAtoms::ATOM, false  );

	// write file
	{
		FILE* fp = fopen( outputGOAFileName, "w" );
		if( fp == 0 ) return false;
		int i;
		for( i=0; i<atomList.size(); i++ )
			printAtomPDB(fp, atomList[i], 0);
		fclose( fp );
	}

	{
		std::vector<PDBParser::Atom*>::iterator iter = atomList.begin(), end = atomList.end();
		for(;iter != end; ++iter)
			delete *iter;
		atomList.clear();
	}
	return true;
}

bool Server::mergeGeometry( int argc, char* argv[] )
{
	if( argc < 5) 
	{
		printUsage();
		return false;
	}

	char inputSurfaceFileName1[256];
	char inputSurfaceFileName2[256];
	char outputSurfaceFileName[256];


	strcpy( inputSurfaceFileName1, argv[2] );
	strcpy( inputSurfaceFileName2, argv[3] );
	strcpy( outputSurfaceFileName, argv[4] );

	Geometry* geometry1 = 0;
	Geometry* geometry2 = 0;

/*	{
		GeometryLoader* gLoader = new GeometryLoader();
		geometry1 = gLoader->loadFile( inputSurfaceFileName1 );
		delete gLoader;
		if( !geometry1 ) return false;
	}

	{
		GeometryLoader* gLoader = new GeometryLoader();
		geometry2 = gLoader->loadFile( inputSurfaceFileName2 );
		delete gLoader;
		if( !geometry2 ) 
		{
			delete geometry1;
			return false;
		}
	}

	Geometry* mergedGeometry = geometry1->merge( geometry2);
	if( !mergedGeometry ) 
	{
		delete geometry1;
		delete geometry2;
		return false;
	}
*/
	Geometry* mergedGeometry = 0;
	{
		GeometryLoader* gLoader = new GeometryLoader();
		mergedGeometry = gLoader->loadFile( "Z:/NMJDatasets/1C2B/1C2B_merged.rawnc" );
		delete gLoader;
		if( !mergedGeometry ) return false;
	}

	Geometry* tg = 0;
	{
		int i;
		for( i=2; i<=140; i++ )
		{
			GeometryLoader* gLoader = new GeometryLoader();
			sprintf( inputSurfaceFileName2, "Z:/NMJDatasets/1C2B/1C2B_%d.rawnc", i );
			geometry1 = gLoader->loadFile( inputSurfaceFileName2 );
			delete gLoader;
			if( !geometry1 ) 
			{
				delete mergedGeometry;
				return false;
			}

			tg = mergedGeometry->merge( geometry1);
			delete mergedGeometry; mergedGeometry = 0;
			delete geometry1; geometry1 = 0;

			mergedGeometry = tg;
			printf("Done %d out of 140 of 1c2b\n", i);
		}
	}

	{
		int i;
		for( i=0; i<=337; i++ )
		{
			GeometryLoader* gLoader = new GeometryLoader();
			sprintf( inputSurfaceFileName2, "Z:/NMJDatasets/2BG9/2BG9_%d.rawnc", i );
			geometry1 = gLoader->loadFile( inputSurfaceFileName2 );
			delete gLoader;
			if( !geometry1 ) 
			{
				delete mergedGeometry;
				return false;
			}

			tg = mergedGeometry->merge( geometry1);
			delete mergedGeometry; mergedGeometry = 0;
			delete geometry1; geometry1 = 0;

			mergedGeometry = tg;

			printf("Done %d out of 337 of 2bg9\n", i);
		}
	}

	{
		int i;
		for( i=0; i<=42; i++ )
		{
			GeometryLoader* gLoader = new GeometryLoader();
			sprintf( inputSurfaceFileName2, "Z:/NMJDatasets/2BG9/2BG9_bottom_%d.rawnc", i );
			geometry1 = gLoader->loadFile( inputSurfaceFileName2 );
			delete gLoader;
			if( !geometry1 ) 
			{
				delete mergedGeometry;
				return false;
			}

			tg = mergedGeometry->merge( geometry1);
			delete mergedGeometry; mergedGeometry = 0;
			delete geometry1; geometry1 = 0;

			mergedGeometry = tg;

			printf("Done %d out of 42 of 2bg9_bottom\n", i);
		}
	}

	{
		GeometryLoader* gLoader = new GeometryLoader();
		bool ret = gLoader->saveFile( outputSurfaceFileName, mergedGeometry );
		delete gLoader;
		delete geometry1;
		delete geometry2;
		delete mergedGeometry;

		if( !ret ) return false;
	}

	return true;
}

bool Server::mergeVolumes( int argc, char* argv[] )
{
	if( argc != 5 && argc != 6 && argc != 11) 
	{
		printUsage();
		return false;
	}

	char inputVolumeFileName1[256];
	char inputVolumeFileName2[256];
	char outputVolumeFileName[256];

	strcpy( inputVolumeFileName1, argv[2] );
	strcpy( inputVolumeFileName2, argv[3] );
	strcpy( outputVolumeFileName, argv[4] );

	bool transform = ( argc == 5 ) ? false : true;

	SimpleVolumeData* sData1 = 0;
	SimpleVolumeData* sData2 = 0;
	
	{
		VolumeLoader* volumeLoader1 = new VolumeLoader();
		sData1 = volumeLoader1->loadFile( inputVolumeFileName1 );
		delete volumeLoader1;
		if( !sData1 ) return false;
		/*//// only for NMJ project! //////
		{
			float minExt[3];
			minExt[0] = sData1->getMinX() * 150;
			minExt[1] = sData1->getMinY() * 150;
			minExt[2] = sData1->getMinZ() * 150;
			sData1->setMinExtent( minExt );
		}
		{
			float maxExt[3];
			maxExt[0] = sData1->getMaxX() * 150;
			maxExt[1] = sData1->getMaxY() * 150;
			maxExt[2] = sData1->getMaxZ() * 150;
			sData1->setMaxExtent( maxExt );
		}
		/////////////////////////////////
*/
		VolumeLoader* volumeLoader2 = new VolumeLoader();
		sData2 = volumeLoader2->loadFile( inputVolumeFileName2 );
		delete volumeLoader2;
		if( !sData2 ) return false;
	}

	FILE* fp = fopen( argv[5], "r" );
	float c1x, c1y, c1z, n1x, n1y, n1z;
	float c2x, c2y, c2z, n2x, n2y, n2z;

	MergeVolumes* mVolumes = new MergeVolumes();

	int numMolecules = 0;
	fscanf(fp, "%d\n", &numMolecules);
	int i;
	for( i=0; i<numMolecules; i++ )
	{
		if( !fscanf( fp, "%f %f %f %f %f %f %f %f %f %f %f %f\n", 
			&c1x, &c1y, &c1z, &n1x, &n1y, &n1z, 
			&c2x, &c2y, &c2z, &n2x, &n2y, &n2z ) ) break;
		mVolumes->setTransformation(c1x, c1y, c1z, n1x, n1y, n1z, 
									c2x, c2y, c2z, n2x, n2y, n2z );

		mVolumes->mergeVolumes( sData1, sData2 );
	}
	
	delete mVolumes;
	fclose( fp );

	{
/*		//// only for NMJ project! //////
		{
			float minExt[3];
			minExt[0] = sData1->getMinX() / 150.0;
			minExt[1] = sData1->getMinY() / 150.0;
			minExt[2] = sData1->getMinZ() / 150.0;
			sData1->setMinExtent( minExt );
		}
		{
			float maxExt[3];
			maxExt[0] = sData1->getMaxX() / 150.0;
			maxExt[1] = sData1->getMaxY() / 150.0;
			maxExt[2] = sData1->getMaxZ() / 150.0;
			sData1->setMaxExtent( maxExt );
		}
*/		/////////////////////////////////
		VolumeLoader* volumeLoader = new VolumeLoader();
		volumeLoader->saveFile(outputVolumeFileName, sData1 );
		delete volumeLoader;
	}
	return true;
}


bool Server::writeTorsionAngles( int argc, char* argv[] )
{
	if( argc != 4 ) 
	{
		printUsage();
		return false;
	}

	PDBParser::GroupOfAtoms* molecule = 0;
	/////// read file
	{
		GOALoader* gLoader = new GOALoader();
		molecule = gLoader->loadFile( argv[2] );
		delete gLoader;
		if( !molecule ) return false;
	}

	FILE* fp = fopen( argv[3], "a");
	if( !fp ) 
	{
		delete molecule;
		return false;
	}

	bool ret = PDBParser::writeTorsionAngles( fp, molecule);

	delete molecule;
	fclose( fp );

	return ret;
}

bool Server::morph( int argc, char* argv[] )
{
	if( argc != 7 ) 
	{
		printUsage();
		return false;
	}

	PDBParser::GroupOfAtoms* molecule1 = 0;
	/////// read file
	{
		GOALoader* gLoader = new GOALoader();
		molecule1 = gLoader->loadFile( argv[2] );
		delete gLoader;
		if( !molecule1 ) return false;
	}

	PDBParser::GroupOfAtoms* molecule2 = 0;
	/////// read file
	{
		GOALoader* gLoader = new GOALoader();
		molecule2 = gLoader->loadFile( argv[3] );
		delete gLoader;
		if( !molecule2 ) { delete molecule1; return false; }
	}
	
	// outputFileNamePrefix = argv[4]
	double resolution = atof( argv[5] );
	int maxSteps = atoi( argv[6] );

	PDBParser::MoleculeMorph* moleculeMorpher = new PDBParser::MoleculeMorph();
	bool ret = moleculeMorpher->morph( molecule1, molecule2, argv[4], resolution, maxSteps );
	delete moleculeMorpher;

	return ret;
}

bool Server::getElecOnSurface( int argc, char* argv[] )
{
	if( argc != 9 ) 
	{
		printUsage();
		return false;
	} 

	VolumeLoader* vLoader = new VolumeLoader();
	SimpleVolumeData* sData = vLoader->loadFile( argv[2] );
	delete vLoader;

	if( !sData ) return false;

	double negCutoff = atof( argv[7] );
	double posCutoff = atof( argv[8] );

	GeometryLoader* gLoader = new GeometryLoader();
	Geometry* inputGeometry = gLoader->loadFile( argv[3] );
	delete gLoader;
	if( !inputGeometry ) return false;
	inputGeometry->CalculateTriSmoothNormals();
	if( !inputGeometry->m_TriVertColors ) inputGeometry->AllocateTriVertColors();

	////// get the function values for each vertex ////////
	double* funcVals = new double[inputGeometry->m_NumTriVerts];
	{
		if( !SimpleVolumeDataIsocontourer::getFunctionValues( sData, inputGeometry, funcVals ) )
		{
			delete []funcVals;
			delete sData;
			delete inputGeometry;
			return false;
		}
	}
	delete sData; sData = 0;
	/////////////////////////////////////////////////////////////


	///// write 3 surfaces depending on the cutoffs ///////////
	Geometry* posGeometry = 0;
	Geometry* neuGeometry = 0;
	Geometry* negGeometry = 0;

	// first get the number of triangles in each
	int numPosTri = 0;
	int numNeuTri = 0;
	int numNegTri = 0;

	{
		int i;
		for( i=0; i<inputGeometry->m_NumTris; i++ )
		{
			double f1 = funcVals[inputGeometry->m_Tris[i*3+0]];
			double f2 = funcVals[inputGeometry->m_Tris[i*3+1]];
			double f3 = funcVals[inputGeometry->m_Tris[i*3+2]];

			if( f1+f2+f3 > 3*posCutoff ) numPosTri++;
			else if( f1+f2+f3 < 3*negCutoff ) numNegTri++;
			else numNeuTri++;
		}
	}

	printf("Number of tris obtained were %d %d %d\n", numPosTri, numNeuTri, numNegTri );

	if( !numPosTri && !numNeuTri && !numNegTri )
	{
		delete []funcVals;
		delete inputGeometry;
		return false;
	}

	if( numPosTri ) posGeometry = new Geometry;
	if( numNeuTri ) neuGeometry = new Geometry;
	if( numNegTri ) negGeometry = new Geometry;

	if( numPosTri ) posGeometry->AllocateTris( inputGeometry->m_NumTriVerts, numPosTri );
	if( numNeuTri ) neuGeometry->AllocateTris( inputGeometry->m_NumTriVerts, numNeuTri );
	if( numNegTri ) negGeometry->AllocateTris( inputGeometry->m_NumTriVerts, numNegTri );

	printf("Allocated all three geometry\n");

	// copy all vertices and relevant triangles to the geometries
	{
		int i;
		for( i=0; i<inputGeometry->m_NumTriVerts*3; i++ )
		{
			if( numPosTri ) posGeometry->m_TriVerts[i] = inputGeometry->m_TriVerts[i];
			if( numNeuTri ) neuGeometry->m_TriVerts[i] = inputGeometry->m_TriVerts[i];
			if( numNegTri ) negGeometry->m_TriVerts[i] = inputGeometry->m_TriVerts[i];
		}
		for( i=0; i<inputGeometry->m_NumTriVerts*3; i++ )
		{
			if( numPosTri ) posGeometry->m_TriVertNormals[i] = inputGeometry->m_TriVertNormals[i];
			if( numNeuTri ) neuGeometry->m_TriVertNormals[i] = inputGeometry->m_TriVertNormals[i];
			if( numNegTri ) negGeometry->m_TriVertNormals[i] = inputGeometry->m_TriVertNormals[i];
		}
		for( i=0; i<inputGeometry->m_NumTriVerts; i++ )
		{
			float r, g, b;
			float f = funcVals[i];

			if( f < negCutoff ) { r = 1.0; g = b = 0; }
			else if( f < posCutoff ) { r = 0.0; g = 0.5; b = 1.0; }
			else if( f < 0 ) { r = 1.0; g = b = 1.0f - f/negCutoff; }
			else if( f > 0 ) { b = 1.0; r = f/posCutoff; g = 0.5 + f / (2.0f*posCutoff); }
			inputGeometry->m_TriVertColors[i*3+0] = r;
			inputGeometry->m_TriVertColors[i*3+1] = g;
			inputGeometry->m_TriVertColors[i*3+2] = b;
		}
		printf("Done with vertices\n");

		int curPosTriIndex = 0;
		int curNeuTriIndex = 0;
		int curNegTriIndex = 0;
		for( i=0; i<inputGeometry->m_NumTris; i++ )
		{
			double f1 = funcVals[inputGeometry->m_Tris[i*3+0]];
			double f2 = funcVals[inputGeometry->m_Tris[i*3+1]];
			double f3 = funcVals[inputGeometry->m_Tris[i*3+2]];

			//printf("Doing with triangle %d with indices [%d %d %d]\n", i, inputGeometry->m_Tris[i*3+0],inputGeometry->m_Tris[i*3+1],inputGeometry->m_Tris[i*3+2] );
			if( f1+f2+f3 > 3*posCutoff ) 
			{
				posGeometry->m_Tris[3*curPosTriIndex+0] = inputGeometry->m_Tris[3*i+0];
				posGeometry->m_Tris[3*curPosTriIndex+1] = inputGeometry->m_Tris[3*i+1];
				posGeometry->m_Tris[3*curPosTriIndex+2] = inputGeometry->m_Tris[3*i+2];
				curPosTriIndex++;
			}
			else if( f1+f2+f3 < 3*negCutoff ) 
			{
				negGeometry->m_Tris[3*curNegTriIndex+0] = inputGeometry->m_Tris[3*i+0];
				negGeometry->m_Tris[3*curNegTriIndex+1] = inputGeometry->m_Tris[3*i+1];
				negGeometry->m_Tris[3*curNegTriIndex+2] = inputGeometry->m_Tris[3*i+2];
				curNegTriIndex++;
			}
			else 
			{
				neuGeometry->m_Tris[3*curNeuTriIndex+0] = inputGeometry->m_Tris[3*i+0];
				neuGeometry->m_Tris[3*curNeuTriIndex+1] = inputGeometry->m_Tris[3*i+1];
				neuGeometry->m_Tris[3*curNeuTriIndex+2] = inputGeometry->m_Tris[3*i+2];
				curNeuTriIndex++;
			}
		}				
	}

	printf("Trying to save\n");
	// save them all
	GeometryLoader* geometryLoader = new GeometryLoader();
	if( numPosTri ) geometryLoader->saveFile( argv[4], posGeometry );
	if( numNeuTri ) geometryLoader->saveFile( argv[5], neuGeometry );
	if( numNegTri ) geometryLoader->saveFile( argv[6], negGeometry );
	geometryLoader->saveFile( argv[3], inputGeometry );
	delete geometryLoader;
	
	if( numPosTri ) delete posGeometry;
	if( numNeuTri ) delete neuGeometry;
	if( numNegTri ) delete negGeometry;
	delete inputGeometry; inputGeometry = 0;

	return true;
}

bool Server::getCurvaturesOnSurface( int argc, char* argv[] )
{
	if( argc != 4 ) 
	{
		printUsage();
		return false;
	} 

	GeometryLoader* gLoader = new GeometryLoader();
	Geometry* inputGeometry = gLoader->loadFile( argv[2] );
	delete gLoader;
	if( !inputGeometry ) return false;
	inputGeometry->computeDerivatives();

	bool ret = inputGeometry->printDerivatives( argv[3] );
	delete inputGeometry;

	return ret;
}
