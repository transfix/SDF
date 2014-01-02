/*****************************************************************************/
/*                                                                           */
/*   Blurmaps, Create volumes, curvatures, surfaces from union of balls      */
/*                                                                           */
/*   Copyright (C) The University of Texas at Austin                         */
/*                                                                           */
/*     Author:     Vinay Siddavanahalli <skvinay@cs.utexas.edu>   2004-2005  */
/*     Author:     John Wiggins         <prok@ices.utexas.edu>    2004-2005  */
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
/*                                                                           */
/*   This library is free software; you can redistribute it and/or           */
/*   modify it under the terms of the GNU Lesser General Public              */
/*   License as published by the Free Software Foundation; either            */
/*   version 2.1 of the License, or (at your option) any later version.      */
/*   Specifically, this library is free for academic or personal non-profit  */
/*   use, with due acknowledgement. Any or all personal profit / industrial  */
/*   use needs to get a proper license approved from us.                     */
/*                                                                           */
/*   This library is distributed in the hope that it will be useful,         */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       */
/*   Lesser General Public License for more details.                         */
/*                                                                           */
/*   You should have received a copy of the GNU Lesser General Public        */
/*   License along with this library; if not, write to the Free Software     */
/*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307    */
/*   USA                                                                     */
/*                                                                           */
/*****************************************************************************/
// SimpleVolumeData.h: interface for the SimpleVolumeData class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CCV_SIMPLE_VOLUME_DATA_H
#define CCV_SIMPLE_VOLUME_DATA_H

#if !WIN32

#ifndef Q_ULLONG
typedef unsigned long long Q_ULLONG;
#endif

#else

#ifndef Q_ULLONG
typedef unsigned __int64 Q_ULLONG;
#endif

#endif

class SimpleVolumeData
{
public:
		SimpleVolumeData(unsigned int dims[3]);
		SimpleVolumeData(unsigned int dim1, unsigned int dim2, unsigned int dim3 );
		~SimpleVolumeData();

		enum DataType {UCHAR=1,USHORT,ULONG,FLOAT,DOUBLE, NO_TYPE=255};

		bool load( const char* fname );
		SimpleVolumeData* createDepthColoredVolume( double* colorMap, int colorMapSize );

		void setNumberOfVariables( unsigned int num );
		void setDimensions( unsigned int dims[3] );
		void setMinExtent( float minExt[3] );
		void setMaxExtent( float maxExt[3] );
		void setData( unsigned int variable, void* data );
		void setType( unsigned int variable, DataType type );
		void setName( unsigned int variable, const char* name );

#ifdef USE_MPI
		void parallelAdjustDimsMerged();
		void parallelAdjustDimsUnmerged();
		unsigned int parallelGetZOffset() const { return m_ZOff; }
		unsigned int parallelGetDepth() const { return m_Dims[3]; }

#endif

		void removeVariable( unsigned int variable );

		bool readRawVFile(const char *filename);
		bool readRawIVFile(const char *filename);

		bool createRawVHeader(char **header, unsigned int *hsize);
		bool createRawIVHeader(char *header);

		void makeVariablesBigEndian();
		void setVariablesToZero();

		unsigned int getWidth() { return m_Dims[0]; }
		unsigned int getHeight() { return m_Dims[1]; }
		unsigned int getDepth() { return m_Dims[2]; }

		float getSpanX() { return (m_Max[0] - m_Min[0]) / (float)(m_Dims[0]-1); }
		float getSpanY() { return (m_Max[1] - m_Min[1]) / (float)(m_Dims[1]-1); }
		float getSpanZ() { return (m_Max[2] - m_Min[2]) / (float)(m_Dims[2]-1); }
		
		float getMinX() { return m_Min[0]; }
		float getMinY() { return m_Min[1]; }
		float getMinZ() { return m_Min[2]; }
		float getMaxX() { return m_Max[0]; }
		float getMaxY() { return m_Max[1]; }
		float getMaxZ() { return m_Max[2]; }

		void* getData( unsigned int variable );
		char** getVariableNames();

		Q_ULLONG getDataSize( unsigned int variable );
		int getTypeSize( unsigned int variable );
		DataType getType( unsigned int variable );
		unsigned int getNumberOfVariables();

		bool getMinMax( unsigned int variable, double* minVal, double* maxVal );
		double getValueAt( unsigned int variable, double x, double y, double z);
		double getValueAt( unsigned int variable, Q_ULLONG index);
		bool setValueAt(unsigned int variable, Q_ULLONG index, double val);
		bool inVolume( double x, double y, double z );

		static bool getMinMax( unsigned char* data, int width, int height, int depth, unsigned char* minVal, unsigned char* maxVal );
		static bool getMinMax( unsigned short* data, int width, int height, int depth, unsigned short* minVal, unsigned short* maxVal );
		static bool getMinMax( unsigned long* data, int width, int height, int depth, unsigned long* minVal, unsigned long* maxVal );
		static bool getMinMax( float* data, int width, int height, int depth, float* minVal, float* maxVal );
		static bool getMinMax( double* data, int width, int height, int depth, double* minVal, double* maxVal );

		static bool getUnsignedCharData( void* data1, int type, unsigned char* data, int width, int height, int depth );
		static bool getNormalizedUnsignedCharData( void* data1, int type, unsigned char* data, int width, int height, int depth );
		static bool getUnsignedCharData( void* data1, int type1, void* data2, int type2, 
										 void* data3, int type3, void* data4, int type4, 
										 unsigned char* data, int width, int height, int depth );
		bool addVolume(SimpleVolumeData* sData, double scale, double sum);
	protected:
		void setDefaults();
		double getDistanceOfVoxel(int i, int j, int k, int width, int height, int depth);

		unsigned char** m_Data;	//	vector valued field

		// general characteristics of data
		DataType* m_DataTypes;
		char** m_VarNames;
		unsigned int m_Dims[4];
		float m_Min[3];
		float m_Max[3];
		unsigned int m_NumberOfVariables;
#ifdef USE_MPI
		unsigned int m_ZOff;
		bool m_MergeOutput;
#endif
};

#endif // CCV_SIMPLE_VOLUME_DATA_H

