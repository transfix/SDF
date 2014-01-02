#include <stdio.h>
#include <string.h>
#include <SimpleVolumeData.h>
#include <sdfLib.h>
#include <Geometry.h>
#include <GeometryLoader.h>
#include <VolumeLoader.h>
#include <cstdlib>

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

int main(int argc, char **argv)
{
  char inputFileName[256];
  char outputVolumeFileName[256];
  int size = 32;
  int flipNormals = 0;

  if(argc != 5 && argc != 11)
   {
     printf("Usage: %s <input geometry> <output volume> <dim (volume size will be dim^3)> <flipNormals (0 or 1)>\n",argv[0]);
     printf("Usage: %s <input geometry> <output volume> <dim (volume size will be dim^3)> <flipNormals (0 or 1)> <minx> <miny> <minz> <maxx> <maxy> <maxz>\n",argv[0]);
     return 0;
   }

  strcpy( inputFileName, argv[1] );
  strcpy( outputVolumeFileName, argv[2] );
  size = atoi(argv[3]);
  flipNormals = atoi(argv[4]);

  char inputFileNameExtension[256];

  if( !getExtension( argv[1], inputFileNameExtension ) )
   {
     printf("Unsupported geometry file.\n");
     return 0;
   }

  GeometryLoader *gLoader = new GeometryLoader();
  Geometry *geometry;
  if( gLoader->isValidExtension( inputFileNameExtension ) )
   {
     GeometryLoader* gLoader = new GeometryLoader();
     geometry = gLoader->loadFile( argv[1] );
     delete gLoader;
   }
 
  if(!geometry)
   {
     printf("Unable to load geometry.\n");
     return 0;
   } 

  float mins[3];
  float maxs[3];

  if(argc == 11) //use extents from command line
    {
      mins[0] = atof(argv[5]); mins[1] = atof(argv[6]); mins[2] = atof(argv[7]);
      maxs[0] = atof(argv[8]); maxs[1] = atof(argv[9]); maxs[2] = atof(argv[10]);
    }
  else
    {
      geometry->CalculateExtents();
      mins[0] = geometry->m_Min[0]; mins[1] = geometry->m_Min[1]; mins[2] = geometry->m_Min[2];
      maxs[0] = geometry->m_Max[0]; maxs[1] = geometry->m_Max[1]; maxs[2] = geometry->m_Max[2];
    }

  SDFLibrary::setParameters(size, flipNormals, mins, maxs);

  //Finally, call the SDF code.
  float* values = SDFLibrary::computeSDF(geometry->m_NumTriVerts, geometry->m_TriVerts, geometry->m_NumTris, (int*)(geometry->m_Tris));
  delete geometry;
  if( !values ) 
   {
     printf("SDFLibrary::computeSDF() failed.\n");
     return 0;
   } 

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
