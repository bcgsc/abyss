#include "HDF5Interface.h"

// dimension names
static const char* X_NAME = "x";
static const char* Y_NAME = "y";
static const char* Z_NAME = "z";
static const char* W_NAME = "w";
static const char* REC_NAME = "sequence";
static const char* SEQ_NAME = "seq";
static const char* STR_NAME = "string";



HDF5Interface::HDF5Interface(const char* filename, FileMode mode) : m_mode(mode), m_recordCounts(SEQLEN, Count3D(SEQLEN, Count2D(SEQLEN, Count1D(SEQLEN))))
{

	// open a netcdf file
	if(m_mode == FM_WRITE)
	{
		openNewFileForWrite(filename);
		printf("file opened\n");
	}
	else
	{
		openFileForRead(filename);
	}
}



void HDF5Interface::openNewFileForWrite(const char* filename)
{
	
	/* Define field information */
	const char *field_names[NFIELDS]  = { "ID","Sequence"};
	hid_t      field_type[NFIELDS];
	hid_t      string_type;
	hsize_t    chunk_size = 10;
	int        *fill_data = NULL;
	int        compress  = 0;
	herr_t     status;
	int        i;
  
  
	/* Initialize field_type */
  	string_type = H5Tcopy( H5T_C_S1 );
	H5Tset_size( string_type, SEQLEN );
	field_type[0] = H5T_NATIVE_INT;
	field_type[1] = string_type;

  
	/* Calculate the size and the offsets of our struct members in memory */
	m_DSTSize =  sizeof( SeqRecord );

	// Calculate offsets 
	m_DSTOffset[0] = HOFFSET( SeqRecord, id );
    m_DSTOffset[1] = HOFFSET( SeqRecord, seq );


	SeqRecord dummy;

	// Calculate field sizes
	m_DSTSizes[0] = sizeof( dummy.id);
	m_DSTSizes[1] = sizeof( dummy.seq);

   	
   	// create file
	m_fileID = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// create tables
	for(int x = 0; x < SEQLEN; x++)
	{
		printf("x=%d\n", x);
		for(int y = 0; y < SEQLEN; y++)
			for(int z = 0; z < SEQLEN; z++)
				for(int w = 0; w < SEQLEN; w++)
				{
					char nameBuffer[64];
					sprintf(nameBuffer, "%d_%d_%d_%d", x,y,z,w);
					
					status = H5TBmake_table( "SeqTable", m_fileID, nameBuffer ,NFIELDS,0,
				                         m_DSTSize,field_names, m_DSTOffset, field_type,
				                         chunk_size, fill_data, compress, NULL);
				}
	}
}

void HDF5Interface::openFileForRead(const char* filename)
{

	
}

void HDF5Interface::writeSequence(Sequence& seq)
{
	const int NUM_RECORDS = 1;
	SeqRecord seqs;
	seqs.id = 0;
	strncpy(seqs.seq,seq.c_str(), SEQLEN);
	
 	int status=H5TBappend_records(m_fileID, "0_0_0_0",NUM_RECORDS, m_DSTSize, m_DSTOffset, m_DSTSizes, &seqs );	
}

void HDF5Interface::readSequencesAtCoordinate(Coord4 pos, Coord4 size)
{

}


HDF5Interface::~HDF5Interface()
{
	H5Fclose( m_fileID );
}
