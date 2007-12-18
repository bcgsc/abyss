# GNU Makefile for the LGAP Assembly suite run make on linux systems, gmake on solaris, etc

# CONDITIONALS
ifeq ($(OSNAME),solaris)
MAKE = gmake
64FLAG = -m64
endif

#########################################################################
# DEFINES

BASE_PATH = ..
COMMON_DIR = Common
DATALAYER_DIR = DataLayer
LGAP_DIR = LGAP
PROOFREAD_DIR = ProofReader
PARTITION_DIR = Partition
TRIMMER_DIR = Trimmer

ALL_DIRS = $(COMMON_DIR) $(DATALAYER_DIR) $(LGAP_DIR) $(PROOFREAD_DIR) $(PARTITION_DIR) $(TRIMMER_DIR)

# Source files
COMMON_SRC = $(wildcard $(COMMON_DIR)/*.cpp)
DATALAYER_SRC = $(wildcard $(DATALAYER_DIR)/*.cpp)
LGAP_SRC = $(wildcard $(LGAP_DIR)/*.cpp)
PROOFREAD_SRC = $(wildcard $(PROOFREAD_DIR)/*.cpp)
PARTITION_SRC = $(wildcard $(PARTITION_DIR)/*.cpp)
TRIMMER_SRC = $(wildcard $(TRIMMER_DIR)/*.cpp)

ALL_SRC = $(COMMON_SRC) $(DATALAYER_SRC) $(LGAP_SRC) $(PROOFREAD_SRC) $(PARTITION_SRC) $(TRIMMER_SRC)

# Object files
COMMON_OBJ = $(patsubst %.cpp,%.o,$(COMMON_SRC))
DATALAYER_OBJ = $(patsubst %.cpp,%.o,$(DATALAYER_SRC))
LGAP_OBJ = $(patsubst %.cpp,%.o,$(LGAP_SRC))
PROOFREAD_OBJ = $(patsubst %.cpp,%.o,$(PROOFREAD_SRC))
PARTITION_OBJ = $(patsubst %.cpp,%.o,$(PARTITION_SRC))
TRIMMER_OBJ = $(patsubst %.cpp,%.o,$(TRIMMER_SRC))

ALL_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(LGAP_OBJ) $(PROOFREAD_OBJ) $(PARTITION_OBJ) $(TRIMMER_OBJ)

LGAP_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(LGAP_OBJ)
PROOFREAD_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(PROOFREAD_OBJ)
PARTITION_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(PARTITION_OBJ)
TRIMMER_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(TRIMMER_OBJ)

# Output dir
BIN_PATH = $(BASE_PATH)/bin

#########################################################################
# FLAGS
CC =		g++

INCLUDE_DIR = $(HOME)/include
LIB_DIR = $(HOME)/lib
OPTIMIZE = -O2 
#DEBUG= -g -D DEBUG -Wall -Wno-sign-compare
#PROFILE = -pg

INCLUDES =	-I$(COMMON_DIR) -I$(DATALAYER_DIR) -I$(INCLUDE_DIR)

CPPFLAGS = $(OPTIMIZE) $(DEBUG) $(PROFILE) $(INCLUDES) 
#LDFLAGS = -Wl,-L$(LIB_DIR),-static,-lnetcdf 
#LDFLAGS = -L$(LIB_DIR) -lnetcdf -lhdf5 -lhdf5_hl -lz -Wl -rpath -L$(LIB_DIR)
LDFLAGS = -Wl,-L$(LIB_DIR),--rpath=$(LIB_DIR),-lnetcdf,-lhdf5,-lhdf5_hl,-lz $(PROFILE)

########################################################################
# Rules

%.o: %.cpp
		$(CC) -c $(CPPFLAGS) $< -o $@ 
		
all:
	$(MAKE) LGAP ProofReader Trimmer Partition 
	
LGAP:		$(LGAP_OBJECTS)
	$(CC) -o $(BIN_PATH)/LGAP $(LGAP_OBJECTS)  $(LDFLAGS)
	
ProofReader:		$(PROOFREAD_OBJECTS)
	$(CC) -o $(BIN_PATH)/ProofReader $(PROOFREAD_OBJECTS) $(LDFLAGS)
	
Partition:		$(PARTITION_OBJECTS)
	$(CC) -o $(BIN_PATH)/Partition $(PARTITION_OBJECTS) $(LDFLAGS)
	
Trimmer:		$(TRIMMER_OBJECTS)
	$(CC) -o $(BIN_PATH)/Trimmer $(TRIMMER_OBJECTS)	$(LDFLAGS)

depend:
	makedepend -Y $(ALL_SRC)

clean:
	$(RM) $(BIN_PATH)/* $(ALL_OBJECTS)


# DO NOT DELETE

Common/HitRecord.o: Common/HitRecord.h Common/Sequence.h Common/CommonDefs.h
Common/HitRecord.o: Common/ReadPrb.h Common/Prb.h
Common/PackedSeq.o: Common/PackedSeq.h Common/Sequence.h Common/CommonDefs.h
Common/PackedSeq.o: Common/ReadPrb.h Common/Prb.h
Common/PairRecord.o: Common/PairRecord.h Common/Sequence.h
Common/PairRecord.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/PhaseSpace.o: Common/PhaseSpace.h Common/Sequence.h
Common/PhaseSpace.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/PhaseSpace.o: Common/HitRecord.h
Common/Prb.o: Common/Prb.h
Common/ReadPrb.o: Common/ReadPrb.h Common/Prb.h
Common/SeqRecord.o: Common/SeqRecord.h Common/CommonDefs.h Common/ReadPrb.h
Common/SeqRecord.o: Common/Prb.h Common/Sequence.h
Common/Sequence.o: Common/Sequence.h Common/CommonDefs.h Common/ReadPrb.h
Common/Sequence.o: Common/Prb.h
Common/SequencePair.o: Common/SequencePair.h Common/Sequence.h
Common/SequencePair.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
DataLayer/HDF5Interface.o: DataLayer/HDF5Interface.h Common/PhaseSpace.h
DataLayer/HDF5Interface.o: Common/Sequence.h Common/CommonDefs.h
DataLayer/HDF5Interface.o: Common/ReadPrb.h Common/Prb.h Common/HitRecord.h
DataLayer/NetCDFInterface.o: DataLayer/NetCDFInterface.h Common/PhaseSpace.h
DataLayer/NetCDFInterface.o: Common/Sequence.h Common/CommonDefs.h
DataLayer/NetCDFInterface.o: Common/ReadPrb.h Common/Prb.h Common/HitRecord.h
DataLayer/PartitionIO.o: DataLayer/PartitionIO.h Common/CommonDefs.h
DataLayer/PartitionIO.o: Common/ReadPrb.h Common/Prb.h Common/Sequence.h
DataLayer/PartitionIO.o: Common/PackedSeq.h
DataLayer/Reader.o: DataLayer/Reader.h Common/CommonDefs.h Common/ReadPrb.h
DataLayer/Reader.o: Common/Prb.h Common/Sequence.h Common/PhaseSpace.h
DataLayer/Reader.o: Common/HitRecord.h Common/PackedSeq.h
DataLayer/Writer.o: DataLayer/Writer.h
LGAP/Path.o: LGAP/Path.h Common/Sequence.h Common/CommonDefs.h
LGAP/Path.o: Common/ReadPrb.h Common/Prb.h Common/SeqRecord.h
LGAP/PathDriver.o: LGAP/PathDriver.h LGAP/Path.h Common/Sequence.h
LGAP/PathDriver.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
LGAP/PathDriver.o: Common/SeqRecord.h Common/HitRecord.h Common/PhaseSpace.h
LGAP/PathDriver.o: Common/PairRecord.h DataLayer/Writer.h
LGAP/PathDriver.o: Common/SequencePair.h
LGAP/PathWalker.o: Common/Sequence.h Common/CommonDefs.h Common/ReadPrb.h
LGAP/PathWalker.o: Common/Prb.h DataLayer/Reader.h Common/PhaseSpace.h
LGAP/PathWalker.o: Common/HitRecord.h Common/PackedSeq.h LGAP/PathDriver.h
LGAP/PathWalker.o: LGAP/Path.h Common/SeqRecord.h Common/PairRecord.h
LGAP/PathWalker.o: DataLayer/Writer.h Common/SequencePair.h
ProofReader/ProofReader.o: ProofReader/ProofReader.h DataLayer/Reader.h
ProofReader/ProofReader.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
ProofReader/ProofReader.o: Common/Sequence.h Common/PhaseSpace.h
ProofReader/ProofReader.o: Common/HitRecord.h Common/PackedSeq.h
ProofReader/ProofReader.o: Common/SeqRecord.h
Partition/Partition.o: DataLayer/Reader.h Common/CommonDefs.h
Partition/Partition.o: Common/ReadPrb.h Common/Prb.h Common/Sequence.h
Partition/Partition.o: Common/PhaseSpace.h Common/HitRecord.h
Partition/Partition.o: Common/PackedSeq.h Partition/Partition.h
Partition/Partition.o: DataLayer/PartitionIO.h
Trimmer/Trimmer.o: Trimmer/Trimmer.h Common/PhaseSpace.h Common/Sequence.h
Trimmer/Trimmer.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Trimmer/Trimmer.o: Common/HitRecord.h DataLayer/Reader.h Common/PackedSeq.h
