#include "NetCDFInterface.h"
#define SEQLEN 36
// dimension names
static const char* X_NAME = "x";
static const char* Y_NAME = "y";
static const char* Z_NAME = "z";
static const char* W_NAME = "w";
static const char* REC_NAME = "sequence";
static const char* SEQ_NAME = "seq";
static const char* STR_NAME = "string";

NetCDFInterface::NetCDFInterface(const char* filename, FileMode mode) : m_mode(mode), m_recordCounts(SEQLEN, Count3D(SEQLEN, Count2D(SEQLEN, Count1D(SEQLEN))))
{
	/*
	// create the record counts structure;
	for(int x = 0; x < SEQLEN; x++)
	{
		m_recordCounts = new int[SEQLEN];
		for(int y = 0; y < SEQLEN; y++)
		{
			m_recordCounts[x] = new
	*/
	
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

void NetCDFInterface::openNewFileForWrite(const char* filename)
{
	int x_varid, y_varid, z_varid, w_varid;
	int x_dimid, y_dimid, z_dimid, w_dimid, str_dimid, rec_dimid;
	int retval;
		
	// open the file for writing and set it up
	// create the file
	if((retval = nc_create(filename, NC_CLOBBER|NC_NETCDF4|NC_SHARE, &m_ncid)))
	{
		ERR(retval);
	}
	
	// define the dimensions
	// the record dimension has unlimited length
	if((retval = nc_def_dim(m_ncid, X_NAME, SEQLEN, &x_dimid)))
		ERR(retval)
	
	if((retval = nc_def_dim(m_ncid, Y_NAME, SEQLEN, &y_dimid)))
		ERR(retval)
	
	if((retval = nc_def_dim(m_ncid, Z_NAME, SEQLEN, &z_dimid)))
		ERR(retval)
	
	if((retval = nc_def_dim(m_ncid, W_NAME, SEQLEN, &w_dimid)))
		ERR(retval)
	
	if((retval = nc_def_dim(m_ncid, STR_NAME, SEQLEN, &str_dimid)))
		ERR(retval)
	
	if((retval = nc_def_dim(m_ncid, REC_NAME, 801, &rec_dimid)))
		ERR(retval)
		
	// define the coordinate variables
	if((retval = nc_def_var(m_ncid, X_NAME, NC_INT, 1, &x_dimid, &x_varid)))
		ERR(retval);

	if((retval = nc_def_var(m_ncid, Y_NAME, NC_INT, 1, &y_dimid, &y_varid)))
		ERR(retval);

	if((retval = nc_def_var(m_ncid, Z_NAME, NC_INT, 1, &z_dimid, &z_varid)))
		ERR(retval);

	if((retval = nc_def_var(m_ncid, W_NAME, NC_INT, 1, &w_dimid, &w_varid)))
		ERR(retval);
		
	// set dimensions
	dimids[0] = rec_dimid;
	dimids[1] = x_dimid;
	dimids[2] = y_dimid;
	dimids[3] = z_dimid;
	dimids[4] = w_dimid;
	dimids[5] = str_dimid;
	
	//define netcdf variables, the actual sequence records, this is saved in the class
	if((retval = nc_def_var(m_ncid, SEQ_NAME, NC_CHAR, NDIMS, dimids, &m_seqVarID)))
		ERR(retval);
	
	// end define mode
	if((retval = nc_enddef(m_ncid)))
		ERR(retval);
		
	// set up the dimension grid
	int seqDataGrid[SEQLEN];
	for(int i = 0; i < SEQLEN; i++)
	{
		seqDataGrid[i] = i;	
	}			
		
	// write coordinate data
	if((retval = nc_put_var_int(m_ncid, x_varid, &seqDataGrid[0])))
		ERR(retval);
	
	if((retval = nc_put_var_int(m_ncid, y_varid, &seqDataGrid[0])))
		ERR(retval);
	
	if((retval = nc_put_var_int(m_ncid, z_varid, &seqDataGrid[0])))
		ERR(retval);
	
	if((retval = nc_put_var_int(m_ncid, w_varid, &seqDataGrid[0])))
		ERR(retval);
	
}

void NetCDFInterface::openFileForRead(const char* filename)
{
	int retval;
	
	// open the file for reading
	if ((retval = nc_open(filename, NC_NOWRITE, &m_ncid)))
		ERR(retval);
		
	// get the number of records so far
	int rec_dimid;
	size_t nRecs;

	if ((retval = nc_inq_dimid(m_ncid, REC_NAME, &rec_dimid)))
		ERR(retval);		
		
	if ((retval = nc_inq_dimlen(m_ncid, rec_dimid, &nRecs)))
		ERR(retval);

	printf("File %s has %d records\n", filename, nRecs);	
	
	// Get the varid for the coordinates
	if ((retval = nc_inq_varid(m_ncid, SEQ_NAME, &m_seqVarID)))
		ERR(retval);
	
}

void NetCDFInterface::writeSequence(Sequence& seq)
{
	
	assert(m_mode == FM_WRITE);
	
	int retval;
	size_t start[NDIMS], count[NDIMS];
	
	// calculate the position of the sequence
	Coord4 pos = PhaseSpace::SequenceToCoord4(seq);
		
	// set up the position to write
	start[1] = pos.x;
	start[2] = pos.y;
	start[3] = pos.z;
	start[4] = pos.w;
	start[5] = 0; // always zero

	// set up the size to write
	count[0] = 1;
	count[1] = 1;
	count[2] = 1;
	count[3] = 1;
	count[4] = 1;
	count[5] = SEQLEN; // always seq len
	
	// write this record in the correct position for the bin
	start[0] = m_recordCounts[pos.x][pos.y][pos.z][pos.w];
	//printf("writing (%d %d %d %d) rec num (%d)", start[1], start[2], start[3], start[4], start[0]);
	if((retval = nc_put_vara_text(m_ncid, m_seqVarID, start, count, seq.c_str())))
		ERR(retval);
		
	//printf("...done\n");
	// increment the record count for this bin
	m_recordCounts[pos.x][pos.y][pos.z][pos.w]++;
}

void NetCDFInterface::readSequencesAtCoordinate(Coord4 pos, Coord4 size)
{
	int retval;
	size_t start[NDIMS], count[NDIMS];

	//read the data

	// start of data slice
	start[1] = pos.x;
	start[2] = pos.y;
	start[3] = pos.z;
	start[4] = pos.w;
	start[5] = 0;	// always start the sequence read from byte 0

	// size of data slice
	count[0] = 1; // num records to read
	count[1] = size.x; // x size
	count[2] = size.y; // y size 
	count[3] = size.z; // z size
	count[4] = size.w; // w size
	count[5] = SEQLEN;
	
	char seqIn[SEQLEN];
	
	// perform the read
	for (int rec = 0; rec < 1; rec++)
	{
		start[0] = rec;
		if ((retval = nc_get_vara_text(m_ncid, m_seqVarID, start, count, seqIn)))
			ERR(retval);

		printf("(%d %d %d %d): %s\n", pos.x, pos.y, pos.z, pos.w, seqIn);
   } 	
}


NetCDFInterface::~NetCDFInterface()
{
	int max = 0;
	for(int x = 0; x < SEQLEN; x++)
	{
		for(int y = 0; y < SEQLEN; y++)
			for(int z = 0; z < SEQLEN; z++)
				for(int w = 0; w < SEQLEN; w++)
				{	
					if(m_recordCounts[x][y][z][w] > max)
					{
						max = m_recordCounts[x][y][z][w];
					}
				}
	}
	
	printf("max record: %d\n", max);
	// close file
	int retval;
	if((retval = nc_close(m_ncid)))
		ERR(retval);
}
