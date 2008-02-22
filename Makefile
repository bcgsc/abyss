# GNU Makefile for the ABYSS Assembly suite run make on linux systems, gmake on solaris, etc

#########################################################################
# DEFINES

BASE_PATH = ..
COMMON_DIR = Common
DATALAYER_DIR = DataLayer
ABYSS_DIR = ABYSS
PROOFREAD_DIR = ProofReader
PARTITION_DIR = Partition
TRIMMER_DIR = Trimmer
REMOVE_MULTIPLICITY_DIR = RemoveMultiplicity
TEST_DIR = Test
PARALLEL_DIR = Parallel

ALL_DIRS = $(COMMON_DIR) $(DATALAYER_DIR) $(ABYSS_DIR) $(PROOFREAD_DIR) $(PARTITION_DIR) $(TRIMMER_DIR) $(REMOVE_MULTIPLICITY_DIR) $(TEST_DIR) $(PARALLEL_DIR)

# Source files
COMMON_SRC = $(wildcard $(COMMON_DIR)/*.cpp)
DATALAYER_SRC = $(wildcard $(DATALAYER_DIR)/*.cpp)
ABYSS_SRC = $(wildcard $(ABYSS_DIR)/*.cpp)
PROOFREAD_SRC = $(wildcard $(PROOFREAD_DIR)/*.cpp)
PARTITION_SRC = $(wildcard $(PARTITION_DIR)/*.cpp)
TRIMMER_SRC = $(wildcard $(TRIMMER_DIR)/*.cpp)
REMOVE_MULTIPLICITY_SRC = $(wildcard $(REMOVE_MULTIPLICITY_DIR)/*.cpp)
VALIDATE_SRC = $(TEST_DIR)/ValidatePartition.cpp
ESTIMATE_CONTIG_SRC = $(TEST_DIR)/EstimateContigs.cpp
PARALLEL_SRC = $(wildcard $(PARALLEL_DIR)/*.cpp)

ALL_SRC = $(COMMON_SRC) $(DATALAYER_SRC) $(ABYSS_SRC) $(PROOFREAD_SRC) $(PARTITION_SRC) $(TRIMMER_SRC) $(REMOVE_MULTIPLICITY_SRC) $(PARALLEL_SRC)

# Object files
COMMON_OBJ = $(patsubst %.cpp,%.o,$(COMMON_SRC))
DATALAYER_OBJ = $(patsubst %.cpp,%.o,$(DATALAYER_SRC))
ABYSS_OBJ = $(patsubst %.cpp,%.o,$(ABYSS_SRC))
PROOFREAD_OBJ = $(patsubst %.cpp,%.o,$(PROOFREAD_SRC))
PARTITION_OBJ = $(patsubst %.cpp,%.o,$(PARTITION_SRC))
TRIMMER_OBJ = $(patsubst %.cpp,%.o,$(TRIMMER_SRC))
REMOVE_MULTIPLICITY_OBJ = $(patsubst %.cpp,%.o,$(REMOVE_MULTIPLICITY_SRC))
VALIDATE_OBJ = $(patsubst %.cpp,%.o,$(VALIDATE_SRC))
ESTIMATE_CONTIG_OBJ = $(patsubst %.cpp,%.o,$(ESTIMATE_CONTIG_SRC))
PARALLEL_OBJ = $(patsubst %.cpp,%.o,$(PARALLEL_SRC))

ALL_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(ABYSS_OBJ) $(PROOFREAD_OBJ) $(PARTITION_OBJ) $(TRIMMER_OBJ) $(ESTIMATE_CONTIG_OBJ) $(VALIDATE_OBJ) $(REMOVE_MULTIPLICITY_OBJ) $(PARALLEL_OBJ)

ABYSS_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(ABYSS_OBJ)
PROOFREAD_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(PROOFREAD_OBJ)
PARTITION_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(PARTITION_OBJ)
TRIMMER_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(TRIMMER_OBJ)
VALIDATE_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(VALIDATE_OBJ)
ESTIMATE_CONTIG_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(ESTIMATE_CONTIG_OBJ)
REMOVE_MULTIPLICITY_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(REMOVE_MULTIPLICITY_OBJ)
PARALLEL_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(PARALLEL_OBJ)

# Output dir
BIN_PATH = $(BASE_PATH)/bin

#########################################################################
# FLAGS
#CC = g++
CC =		mpic++

INCLUDE_DIR = $(HOME)/include
LIB_DIR = $(HOME)/lib
OPTIMIZE = -O2
#64BIT = -m64
#DEBUG= -g -D DEBUG -Wall -Wno-sign-compare
#PROFILE = -pg

INCLUDES =	-I$(COMMON_DIR) -I$(DATALAYER_DIR) -I$(INCLUDE_DIR)

CPPFLAGS = $(OPTIMIZE) $(DEBUG) $(PROFILE) $(INCLUDES) $(64BIT)
#LDFLAGS = -Wl,-L$(LIB_DIR),-static,-lnetcdf 
#LDFLAGS = -L$(LIB_DIR) -lnetcdf -lhdf5 -lhdf5_hl -lz -Wl -rpath -L$(LIB_DIR)
LDFLAGS = -L$(LIB_DIR) $(PROFILE) $(64BIT)

########################################################################
# Rules

%.o: %.cpp
		$(CC) -c $(CPPFLAGS) $< -o $@ 
		
all:
	$(MAKE) ABYSS ProofReader Partition Trimmer RemoveMultiplicity test ABYSS-P
	
test:
	$(MAKE) ValidatePartition EstimateContigs
	
ABYSS:		$(ABYSS_OBJECTS)
	$(CC) -o $(BIN_PATH)/ABYSS $(ABYSS_OBJECTS)  $(LDFLAGS)
	
ProofReader:		$(PROOFREAD_OBJECTS)
	$(CC) -o $(BIN_PATH)/ProofReader $(PROOFREAD_OBJECTS) $(LDFLAGS)
	
Partition:		$(PARTITION_OBJECTS)
	$(CC) -o $(BIN_PATH)/Partition $(PARTITION_OBJECTS) $(LDFLAGS)
	
Trimmer:		$(TRIMMER_OBJECTS)
	$(CC) -o $(BIN_PATH)/Trimmer $(TRIMMER_OBJECTS)	$(LDFLAGS)
	
ValidatePartition:		$(VALIDATE_OBJECTS)
	$(CC) -o $(BIN_PATH)/ValidatePartition $(VALIDATE_OBJECTS)	$(LDFLAGS)
	
RemoveMultiplicity:		$(REMOVE_MULTIPLICITY_OBJECTS)
	$(CC) -o $(BIN_PATH)/RemoveMultiplicity $(REMOVE_MULTIPLICITY_OBJECTS)	$(LDFLAGS)
	
EstimateContigs:	$(ESTIMATE_CONTIG_OBJECTS)
	$(CC) -o $(BIN_PATH)/EstimateContigs $(ESTIMATE_CONTIG_OBJECTS)	$(LDFLAGS)
	
ABYSS-P:	$(PARALLEL_OBJECTS)
	$(CC) -o $(BIN_PATH)/ABYSS-P $(PARALLEL_OBJECTS) $(LDFLAGS)

depend:
	makedepend -Y $(ALL_SRC)

clean:
	$(RM) $(BIN_PATH)/* $(ALL_OBJECTS)


# DO NOT DELETE

Common/AssemblyData.o: Common/AssemblyData.h Common/Sequence.h
Common/AssemblyData.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/AssemblyData.o: Common/PackedSeq.h Common/CommonUtils.h
Common/AssemblyData.o: Common/SeqExt.h Common/HitRecord.h
Common/CommonUtils.o: Common/CommonUtils.h Common/PackedSeq.h
Common/CommonUtils.o: Common/Sequence.h Common/CommonDefs.h Common/ReadPrb.h
Common/CommonUtils.o: Common/Prb.h Common/SeqExt.h
Common/Config.o: Common/Config.h
Common/HitRecord.o: Common/HitRecord.h Common/CommonDefs.h Common/ReadPrb.h
Common/HitRecord.o: Common/Prb.h Common/PackedSeq.h Common/CommonUtils.h
Common/HitRecord.o: Common/Sequence.h Common/SeqExt.h
Common/PackedSeq.o: Common/PackedSeq.h Common/CommonUtils.h Common/Sequence.h
Common/PackedSeq.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/PackedSeq.o: Common/SeqExt.h
Common/PairRecord.o: Common/PairRecord.h Common/PackedSeq.h
Common/PairRecord.o: Common/CommonUtils.h Common/Sequence.h
Common/PairRecord.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/PairRecord.o: Common/SeqExt.h
Common/PartitionLoader.o: Common/PartitionLoader.h Common/CommonDefs.h
Common/PartitionLoader.o: Common/ReadPrb.h Common/Prb.h Common/PhaseSpace.h
Common/PartitionLoader.o: Common/Sequence.h Common/PackedSeq.h
Common/PartitionLoader.o: Common/CommonUtils.h Common/SeqExt.h
Common/PartitionLoader.o: Common/HitRecord.h Common/Config.h
Common/PhaseSpace.o: Common/PhaseSpace.h Common/Sequence.h
Common/PhaseSpace.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/PhaseSpace.o: Common/PackedSeq.h Common/CommonUtils.h Common/SeqExt.h
Common/PhaseSpace.o: Common/HitRecord.h
Common/Prb.o: Common/Prb.h
Common/ReadPrb.o: Common/ReadPrb.h Common/Prb.h
Common/SeqExt.o: Common/SeqExt.h Common/CommonDefs.h Common/ReadPrb.h
Common/SeqExt.o: Common/Prb.h
Common/SeqRecord.o: Common/SeqRecord.h Common/PackedSeq.h
Common/SeqRecord.o: Common/CommonUtils.h Common/Sequence.h
Common/SeqRecord.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/SeqRecord.o: Common/SeqExt.h
Common/Sequence.o: Common/Sequence.h Common/CommonDefs.h Common/ReadPrb.h
Common/Sequence.o: Common/Prb.h
Common/SequencePair.o: Common/SequencePair.h Common/PackedSeq.h
Common/SequencePair.o: Common/CommonUtils.h Common/Sequence.h
Common/SequencePair.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/SequencePair.o: Common/SeqExt.h
DataLayer/FastaReader.o: DataLayer/FastaReader.h DataLayer/IFileReader.h
DataLayer/FastaReader.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
DataLayer/FastaReader.o: Common/Sequence.h Common/PackedSeq.h
DataLayer/FastaReader.o: Common/CommonUtils.h Common/SeqExt.h
DataLayer/FastaWriter.o: DataLayer/FastaWriter.h Common/CommonDefs.h
DataLayer/FastaWriter.o: Common/ReadPrb.h Common/Prb.h Common/Sequence.h
DataLayer/FastaWriter.o: DataLayer/IFileWriter.h Common/PackedSeq.h
DataLayer/FastaWriter.o: Common/CommonUtils.h Common/SeqExt.h
DataLayer/PackedSeqReader.o: DataLayer/PackedSeqReader.h
DataLayer/PackedSeqReader.o: DataLayer/IFileReader.h Common/CommonDefs.h
DataLayer/PackedSeqReader.o: Common/ReadPrb.h Common/Prb.h Common/Sequence.h
DataLayer/PackedSeqReader.o: Common/PackedSeq.h Common/CommonUtils.h
DataLayer/PackedSeqReader.o: Common/SeqExt.h
DataLayer/PackedSeqWriter.o: DataLayer/PackedSeqWriter.h Common/CommonDefs.h
DataLayer/PackedSeqWriter.o: Common/ReadPrb.h Common/Prb.h Common/Sequence.h
DataLayer/PackedSeqWriter.o: DataLayer/IFileWriter.h Common/PackedSeq.h
DataLayer/PackedSeqWriter.o: Common/CommonUtils.h Common/SeqExt.h
DataLayer/Reader.o: DataLayer/Reader.h Common/CommonDefs.h Common/ReadPrb.h
DataLayer/Reader.o: Common/Prb.h Common/Sequence.h Common/PhaseSpace.h
DataLayer/Reader.o: Common/PackedSeq.h Common/CommonUtils.h Common/SeqExt.h
DataLayer/Reader.o: Common/HitRecord.h
DataLayer/SequenceCollection.o: DataLayer/SequenceCollection.h
DataLayer/SequenceCollection.o: DataLayer/ISequenceCollection.h
DataLayer/SequenceCollection.o: Common/CommonDefs.h Common/ReadPrb.h
DataLayer/SequenceCollection.o: Common/Prb.h Common/CommonUtils.h
DataLayer/SequenceCollection.o: Common/PackedSeq.h Common/Sequence.h
DataLayer/SequenceCollection.o: Common/SeqExt.h Common/HitRecord.h
DataLayer/Writer.o: DataLayer/Writer.h
ABYSS/Abyss.o: ABYSS/Abyss.h Common/Sequence.h Common/CommonDefs.h
ABYSS/Abyss.o: Common/ReadPrb.h Common/Prb.h DataLayer/Reader.h
ABYSS/Abyss.o: Common/PhaseSpace.h Common/PackedSeq.h Common/CommonUtils.h
ABYSS/Abyss.o: Common/SeqExt.h Common/HitRecord.h ABYSS/PathDriver.h
ABYSS/Abyss.o: ABYSS/Path.h Common/SeqRecord.h Common/PairRecord.h
ABYSS/Abyss.o: DataLayer/Writer.h Common/SequencePair.h Common/Config.h
ABYSS/Abyss.o: Common/PartitionLoader.h DataLayer/FastaWriter.h
ABYSS/Abyss.o: DataLayer/IFileWriter.h DataLayer/FastaReader.h
ABYSS/Abyss.o: DataLayer/IFileReader.h DataLayer/PackedSeqWriter.h
ABYSS/Path.o: ABYSS/Path.h Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
ABYSS/Path.o: Common/SeqRecord.h Common/PackedSeq.h Common/CommonUtils.h
ABYSS/Path.o: Common/Sequence.h Common/SeqExt.h
ABYSS/PathDriver.o: ABYSS/PathDriver.h Common/CommonDefs.h Common/ReadPrb.h
ABYSS/PathDriver.o: Common/Prb.h ABYSS/Path.h Common/SeqRecord.h
ABYSS/PathDriver.o: Common/PackedSeq.h Common/CommonUtils.h Common/Sequence.h
ABYSS/PathDriver.o: Common/SeqExt.h Common/HitRecord.h Common/PhaseSpace.h
ABYSS/PathDriver.o: Common/PairRecord.h DataLayer/Writer.h
ABYSS/PathDriver.o: Common/SequencePair.h
ProofReader/ProofReader.o: ProofReader/ProofReader.h DataLayer/Reader.h
ProofReader/ProofReader.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
ProofReader/ProofReader.o: Common/Sequence.h Common/PhaseSpace.h
ProofReader/ProofReader.o: Common/PackedSeq.h Common/CommonUtils.h
ProofReader/ProofReader.o: Common/SeqExt.h Common/HitRecord.h
ProofReader/ProofReader.o: Common/SeqRecord.h
Partition/Partition.o: DataLayer/Reader.h Common/CommonDefs.h
Partition/Partition.o: Common/ReadPrb.h Common/Prb.h Common/Sequence.h
Partition/Partition.o: Common/PhaseSpace.h Common/PackedSeq.h
Partition/Partition.o: Common/CommonUtils.h Common/SeqExt.h
Partition/Partition.o: Common/HitRecord.h Partition/Partition.h
Partition/Partition.o: DataLayer/IFileWriter.h DataLayer/PackedSeqWriter.h
Partition/Partition.o: Common/Config.h DataLayer/IFileReader.h
Partition/Partition.o: DataLayer/FastaReader.h DataLayer/PackedSeqReader.h
Partition/Partition.o: DataLayer/FastaWriter.h
Trimmer/Trimmer.o: Trimmer/Trimmer.h DataLayer/Reader.h Common/CommonDefs.h
Trimmer/Trimmer.o: Common/ReadPrb.h Common/Prb.h Common/Sequence.h
Trimmer/Trimmer.o: Common/PhaseSpace.h Common/PackedSeq.h
Trimmer/Trimmer.o: Common/CommonUtils.h Common/SeqExt.h Common/HitRecord.h
Trimmer/Trimmer.o: Common/PartitionLoader.h Common/Config.h
Trimmer/Trimmer.o: DataLayer/PackedSeqWriter.h DataLayer/IFileWriter.h
RemoveMultiplicity/RemoveMultiplicity.o: RemoveMultiplicity/RemoveMultiplicity.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/Config.h Common/CommonDefs.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/ReadPrb.h Common/Prb.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/PhaseSpace.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/Sequence.h Common/PackedSeq.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/CommonUtils.h Common/SeqExt.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/HitRecord.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/Reader.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/PartitionLoader.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/PackedSeqWriter.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/IFileWriter.h
Parallel/CommLayer.o: Parallel/CommLayer.h Common/PackedSeq.h
Parallel/CommLayer.o: Common/CommonUtils.h Common/Sequence.h
Parallel/CommLayer.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Parallel/CommLayer.o: Common/SeqExt.h
Parallel/DistributedPhaseSpace.o: Parallel/DistributedPhaseSpace.h
Parallel/DistributedPhaseSpace.o: Parallel/CommLayer.h Common/PackedSeq.h
Parallel/DistributedPhaseSpace.o: Common/CommonUtils.h Common/Sequence.h
Parallel/DistributedPhaseSpace.o: Common/CommonDefs.h Common/ReadPrb.h
Parallel/DistributedPhaseSpace.o: Common/Prb.h Common/SeqExt.h
Parallel/NetworkSequenceSpace.o: Parallel/NetworkSequenceSpace.h
Parallel/NetworkSequenceSpace.o: Common/PackedSeq.h Common/CommonUtils.h
Parallel/NetworkSequenceSpace.o: Common/Sequence.h Common/CommonDefs.h
Parallel/NetworkSequenceSpace.o: Common/ReadPrb.h Common/Prb.h
Parallel/NetworkSequenceSpace.o: Common/SeqExt.h Parallel/CommLayer.h
Parallel/parallelAbyss.o: Parallel/parallelAbyss.h Common/PackedSeq.h
Parallel/parallelAbyss.o: Common/CommonUtils.h Common/Sequence.h
Parallel/parallelAbyss.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Parallel/parallelAbyss.o: Common/SeqExt.h DataLayer/FastaReader.h
Parallel/parallelAbyss.o: DataLayer/IFileReader.h
Parallel/parallelAbyss.o: Parallel/DistributedPhaseSpace.h
Parallel/parallelAbyss.o: Parallel/CommLayer.h
