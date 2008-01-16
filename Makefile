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
REMOVE_MULTIPLICITY_DIR = RemoveMultiplicity
TEST_DIR = Test

ALL_DIRS = $(COMMON_DIR) $(DATALAYER_DIR) $(LGAP_DIR) $(PROOFREAD_DIR) $(PARTITION_DIR) $(TRIMMER_DIR) $(REMOVE_MULTIPLICITY_DIR) $(TEST_DIR)

# Source files
COMMON_SRC = $(wildcard $(COMMON_DIR)/*.cpp)
DATALAYER_SRC = $(wildcard $(DATALAYER_DIR)/*.cpp)
LGAP_SRC = $(wildcard $(LGAP_DIR)/*.cpp)
PROOFREAD_SRC = $(wildcard $(PROOFREAD_DIR)/*.cpp)
PARTITION_SRC = $(wildcard $(PARTITION_DIR)/*.cpp)
TRIMMER_SRC = $(wildcard $(TRIMMER_DIR)/*.cpp)
REMOVE_MULTIPLICITY_SRC = $(wildcard $(REMOVE_MULTIPLICITY_DIR)/*.cpp)
VALIDATE_SRC = $(TEST_DIR)/ValidatePartition.cpp
ESTIMATE_CONTIG_SRC = $(TEST_DIR)/EstimateContigs.cpp

ALL_SRC = $(COMMON_SRC) $(DATALAYER_SRC) $(LGAP_SRC) $(PROOFREAD_SRC) $(PARTITION_SRC) $(TRIMMER_SRC) $(REMOVE_MULTIPLICITY_SRC)

# Object files
COMMON_OBJ = $(patsubst %.cpp,%.o,$(COMMON_SRC))
DATALAYER_OBJ = $(patsubst %.cpp,%.o,$(DATALAYER_SRC))
LGAP_OBJ = $(patsubst %.cpp,%.o,$(LGAP_SRC))
PROOFREAD_OBJ = $(patsubst %.cpp,%.o,$(PROOFREAD_SRC))
PARTITION_OBJ = $(patsubst %.cpp,%.o,$(PARTITION_SRC))
TRIMMER_OBJ = $(patsubst %.cpp,%.o,$(TRIMMER_SRC))
REMOVE_MULTIPLICITY_OBJ = $(patsubst %.cpp,%.o,$(REMOVE_MULTIPLICITY_SRC))
VALIDATE_OBJ = $(patsubst %.cpp,%.o,$(VALIDATE_SRC))
ESTIMATE_CONTIG_OBJ = $(patsubst %.cpp,%.o,$(ESTIMATE_CONTIG_SRC))

ALL_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(LGAP_OBJ) $(PROOFREAD_OBJ) $(PARTITION_OBJ) $(TRIMMER_OBJ) $(ESTIMATE_CONTIG_OBJ) $(VALIDATE_OBJ) $(REMOVE_MULTIPLICITY_OBJ)

LGAP_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(LGAP_OBJ)
PROOFREAD_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(PROOFREAD_OBJ)
PARTITION_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(PARTITION_OBJ)
TRIMMER_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(TRIMMER_OBJ)
VALIDATE_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(VALIDATE_OBJ)
ESTIMATE_CONTIG_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(ESTIMATE_CONTIG_OBJ)
REMOVE_MULTIPLICITY_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(REMOVE_MULTIPLICITY_OBJ)

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
LDFLAGS = -L$(LIB_DIR)

########################################################################
# Rules

%.o: %.cpp
		$(CC) -c $(CPPFLAGS) $< -o $@ 
		
all:
	$(MAKE) LGAP ProofReader Partition Trimmer RemoveMultiplicity test
	
test:
	$(MAKE) ValidatePartition EstimateContigs
	
LGAP:		$(LGAP_OBJECTS)
	$(CC) -o $(BIN_PATH)/LGAP $(LGAP_OBJECTS)  $(LDFLAGS)
	
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

depend:
	makedepend $(ALL_SRC)

clean:
	$(RM) $(BIN_PATH)/* $(ALL_OBJECTS)


# DO NOT DELETE

Common/CommonUtils.o: /usr/include/math.h /usr/include/features.h
Common/CommonUtils.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/CommonUtils.o: /usr/include/bits/huge_val.h
Common/CommonUtils.o: /usr/include/bits/mathdef.h
Common/CommonUtils.o: /usr/include/bits/mathcalls.h Common/CommonUtils.h
Common/CommonUtils.o: Common/PackedSeq.h Common/Sequence.h
Common/CommonUtils.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/CommonUtils.o: /usr/include/assert.h
Common/Config.o: Common/Config.h
Common/HitRecord.o: Common/HitRecord.h Common/CommonDefs.h Common/ReadPrb.h
Common/HitRecord.o: Common/Prb.h /usr/include/assert.h
Common/HitRecord.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/HitRecord.o: /usr/include/gnu/stubs.h Common/PackedSeq.h
Common/HitRecord.o: Common/Sequence.h
Common/PackedSeq.o: /usr/include/math.h /usr/include/features.h
Common/PackedSeq.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/PackedSeq.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
Common/PackedSeq.o: /usr/include/bits/mathcalls.h Common/PackedSeq.h
Common/PackedSeq.o: Common/Sequence.h Common/CommonDefs.h Common/ReadPrb.h
Common/PackedSeq.o: Common/Prb.h /usr/include/assert.h
Common/PairRecord.o: Common/PairRecord.h Common/PackedSeq.h Common/Sequence.h
Common/PairRecord.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/PairRecord.o: /usr/include/assert.h /usr/include/features.h
Common/PairRecord.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/PartitionLoader.o: Common/PartitionLoader.h Common/CommonDefs.h
Common/PartitionLoader.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/PartitionLoader.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/PartitionLoader.o: /usr/include/gnu/stubs.h Common/PhaseSpace.h
Common/PartitionLoader.o: /usr/include/stdio.h
Common/PartitionLoader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
Common/PartitionLoader.o: /usr/include/bits/types.h
Common/PartitionLoader.o: /usr/include/bits/wordsize.h
Common/PartitionLoader.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Common/PartitionLoader.o: /usr/include/_G_config.h /usr/include/wchar.h
Common/PartitionLoader.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Common/PartitionLoader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
Common/PartitionLoader.o: /usr/include/bits/stdio_lim.h
Common/PartitionLoader.o: /usr/include/bits/sys_errlist.h Common/Sequence.h
Common/PartitionLoader.o: Common/PackedSeq.h Common/HitRecord.h
Common/PartitionLoader.o: Common/Config.h
Common/PhaseSpace.o: Common/PhaseSpace.h /usr/include/stdio.h
Common/PhaseSpace.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/PhaseSpace.o: /usr/include/gnu/stubs.h
Common/PhaseSpace.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
Common/PhaseSpace.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
Common/PhaseSpace.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Common/PhaseSpace.o: /usr/include/_G_config.h /usr/include/wchar.h
Common/PhaseSpace.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Common/PhaseSpace.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
Common/PhaseSpace.o: /usr/include/bits/stdio_lim.h
Common/PhaseSpace.o: /usr/include/bits/sys_errlist.h Common/Sequence.h
Common/PhaseSpace.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/PhaseSpace.o: /usr/include/assert.h Common/PackedSeq.h
Common/PhaseSpace.o: Common/HitRecord.h Common/CommonUtils.h
Common/Prb.o: /usr/include/math.h /usr/include/features.h
Common/Prb.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/Prb.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
Common/Prb.o: /usr/include/bits/mathcalls.h Common/Prb.h
Common/Prb.o: /usr/include/assert.h
Common/ReadPrb.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/ReadPrb.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/ReadPrb.o: /usr/include/gnu/stubs.h
Common/SeqRecord.o: Common/SeqRecord.h Common/PackedSeq.h Common/Sequence.h
Common/SeqRecord.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/SeqRecord.o: /usr/include/assert.h /usr/include/features.h
Common/SeqRecord.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/Sequence.o: /usr/include/stdio.h /usr/include/features.h
Common/Sequence.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/Sequence.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
Common/Sequence.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
Common/Sequence.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Common/Sequence.o: /usr/include/_G_config.h /usr/include/wchar.h
Common/Sequence.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Common/Sequence.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
Common/Sequence.o: /usr/include/bits/stdio_lim.h
Common/Sequence.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
Common/Sequence.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
Common/Sequence.o: /usr/include/bits/mathcalls.h Common/Sequence.h
Common/Sequence.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Common/Sequence.o: /usr/include/assert.h
Common/SequencePair.o: Common/SequencePair.h Common/PackedSeq.h
Common/SequencePair.o: Common/Sequence.h Common/CommonDefs.h Common/ReadPrb.h
Common/SequencePair.o: Common/Prb.h /usr/include/assert.h
Common/SequencePair.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/SequencePair.o: /usr/include/gnu/stubs.h
DataLayer/FastaReader.o: DataLayer/FastaReader.h /usr/include/stdio.h
DataLayer/FastaReader.o: /usr/include/features.h /usr/include/sys/cdefs.h
DataLayer/FastaReader.o: /usr/include/gnu/stubs.h
DataLayer/FastaReader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
DataLayer/FastaReader.o: /usr/include/bits/types.h
DataLayer/FastaReader.o: /usr/include/bits/wordsize.h
DataLayer/FastaReader.o: /usr/include/bits/typesizes.h /usr/include/libio.h
DataLayer/FastaReader.o: /usr/include/_G_config.h /usr/include/wchar.h
DataLayer/FastaReader.o: /usr/include/bits/wchar.h /usr/include/gconv.h
DataLayer/FastaReader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
DataLayer/FastaReader.o: /usr/include/bits/stdio_lim.h
DataLayer/FastaReader.o: /usr/include/bits/sys_errlist.h
DataLayer/FastaReader.o: DataLayer/IFileReader.h Common/CommonDefs.h
DataLayer/FastaReader.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
DataLayer/FastaReader.o: Common/Sequence.h Common/PackedSeq.h
DataLayer/FastaWriter.o: DataLayer/FastaWriter.h /usr/include/stdio.h
DataLayer/FastaWriter.o: /usr/include/features.h /usr/include/sys/cdefs.h
DataLayer/FastaWriter.o: /usr/include/gnu/stubs.h
DataLayer/FastaWriter.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
DataLayer/FastaWriter.o: /usr/include/bits/types.h
DataLayer/FastaWriter.o: /usr/include/bits/wordsize.h
DataLayer/FastaWriter.o: /usr/include/bits/typesizes.h /usr/include/libio.h
DataLayer/FastaWriter.o: /usr/include/_G_config.h /usr/include/wchar.h
DataLayer/FastaWriter.o: /usr/include/bits/wchar.h /usr/include/gconv.h
DataLayer/FastaWriter.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
DataLayer/FastaWriter.o: /usr/include/bits/stdio_lim.h
DataLayer/FastaWriter.o: /usr/include/bits/sys_errlist.h Common/CommonDefs.h
DataLayer/FastaWriter.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
DataLayer/FastaWriter.o: Common/Sequence.h DataLayer/IFileWriter.h
DataLayer/FastaWriter.o: Common/PackedSeq.h
DataLayer/PackedSeqReader.o: DataLayer/PackedSeqReader.h /usr/include/stdio.h
DataLayer/PackedSeqReader.o: /usr/include/features.h /usr/include/sys/cdefs.h
DataLayer/PackedSeqReader.o: /usr/include/gnu/stubs.h
DataLayer/PackedSeqReader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
DataLayer/PackedSeqReader.o: /usr/include/bits/types.h
DataLayer/PackedSeqReader.o: /usr/include/bits/wordsize.h
DataLayer/PackedSeqReader.o: /usr/include/bits/typesizes.h
DataLayer/PackedSeqReader.o: /usr/include/libio.h /usr/include/_G_config.h
DataLayer/PackedSeqReader.o: /usr/include/wchar.h /usr/include/bits/wchar.h
DataLayer/PackedSeqReader.o: /usr/include/gconv.h
DataLayer/PackedSeqReader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
DataLayer/PackedSeqReader.o: /usr/include/bits/stdio_lim.h
DataLayer/PackedSeqReader.o: /usr/include/bits/sys_errlist.h
DataLayer/PackedSeqReader.o: DataLayer/IFileReader.h Common/CommonDefs.h
DataLayer/PackedSeqReader.o: Common/ReadPrb.h Common/Prb.h
DataLayer/PackedSeqReader.o: /usr/include/assert.h Common/Sequence.h
DataLayer/PackedSeqReader.o: Common/PackedSeq.h
DataLayer/PackedSeqWriter.o: DataLayer/PackedSeqWriter.h /usr/include/stdio.h
DataLayer/PackedSeqWriter.o: /usr/include/features.h /usr/include/sys/cdefs.h
DataLayer/PackedSeqWriter.o: /usr/include/gnu/stubs.h
DataLayer/PackedSeqWriter.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/types.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/wordsize.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/typesizes.h
DataLayer/PackedSeqWriter.o: /usr/include/libio.h /usr/include/_G_config.h
DataLayer/PackedSeqWriter.o: /usr/include/wchar.h /usr/include/bits/wchar.h
DataLayer/PackedSeqWriter.o: /usr/include/gconv.h
DataLayer/PackedSeqWriter.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/stdio_lim.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/sys_errlist.h
DataLayer/PackedSeqWriter.o: Common/CommonDefs.h Common/ReadPrb.h
DataLayer/PackedSeqWriter.o: Common/Prb.h /usr/include/assert.h
DataLayer/PackedSeqWriter.o: Common/Sequence.h DataLayer/IFileWriter.h
DataLayer/PackedSeqWriter.o: Common/PackedSeq.h
DataLayer/Reader.o: /usr/include/stdio.h /usr/include/features.h
DataLayer/Reader.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
DataLayer/Reader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
DataLayer/Reader.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
DataLayer/Reader.o: /usr/include/bits/typesizes.h /usr/include/libio.h
DataLayer/Reader.o: /usr/include/_G_config.h /usr/include/wchar.h
DataLayer/Reader.o: /usr/include/bits/wchar.h /usr/include/gconv.h
DataLayer/Reader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
DataLayer/Reader.o: /usr/include/bits/stdio_lim.h
DataLayer/Reader.o: /usr/include/bits/sys_errlist.h DataLayer/Reader.h
DataLayer/Reader.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
DataLayer/Reader.o: /usr/include/assert.h Common/Sequence.h
DataLayer/Reader.o: Common/PhaseSpace.h Common/PackedSeq.h Common/HitRecord.h
DataLayer/Writer.o: DataLayer/Writer.h
LGAP/Path.o: LGAP/Path.h Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
LGAP/Path.o: /usr/include/assert.h /usr/include/features.h
LGAP/Path.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
LGAP/Path.o: Common/SeqRecord.h Common/PackedSeq.h Common/Sequence.h
LGAP/Path.o: Common/CommonUtils.h
LGAP/PathDriver.o: LGAP/PathDriver.h Common/CommonDefs.h Common/ReadPrb.h
LGAP/PathDriver.o: Common/Prb.h /usr/include/assert.h /usr/include/features.h
LGAP/PathDriver.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
LGAP/PathDriver.o: LGAP/Path.h Common/SeqRecord.h Common/PackedSeq.h
LGAP/PathDriver.o: Common/Sequence.h Common/HitRecord.h Common/PhaseSpace.h
LGAP/PathDriver.o: /usr/include/stdio.h
LGAP/PathDriver.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
LGAP/PathDriver.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
LGAP/PathDriver.o: /usr/include/bits/typesizes.h /usr/include/libio.h
LGAP/PathDriver.o: /usr/include/_G_config.h /usr/include/wchar.h
LGAP/PathDriver.o: /usr/include/bits/wchar.h /usr/include/gconv.h
LGAP/PathDriver.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
LGAP/PathDriver.o: /usr/include/bits/stdio_lim.h
LGAP/PathDriver.o: /usr/include/bits/sys_errlist.h Common/PairRecord.h
LGAP/PathDriver.o: DataLayer/Writer.h Common/SequencePair.h
LGAP/PathWalker.o: /usr/include/stdio.h /usr/include/features.h
LGAP/PathWalker.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
LGAP/PathWalker.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
LGAP/PathWalker.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
LGAP/PathWalker.o: /usr/include/bits/typesizes.h /usr/include/libio.h
LGAP/PathWalker.o: /usr/include/_G_config.h /usr/include/wchar.h
LGAP/PathWalker.o: /usr/include/bits/wchar.h /usr/include/gconv.h
LGAP/PathWalker.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
LGAP/PathWalker.o: /usr/include/bits/stdio_lim.h
LGAP/PathWalker.o: /usr/include/bits/sys_errlist.h Common/Sequence.h
LGAP/PathWalker.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
LGAP/PathWalker.o: /usr/include/assert.h DataLayer/Reader.h
LGAP/PathWalker.o: Common/PhaseSpace.h Common/PackedSeq.h Common/HitRecord.h
LGAP/PathWalker.o: LGAP/PathDriver.h LGAP/Path.h Common/SeqRecord.h
LGAP/PathWalker.o: Common/PairRecord.h DataLayer/Writer.h
LGAP/PathWalker.o: Common/SequencePair.h
ProofReader/ProofReader.o: ProofReader/ProofReader.h DataLayer/Reader.h
ProofReader/ProofReader.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
ProofReader/ProofReader.o: /usr/include/assert.h /usr/include/features.h
ProofReader/ProofReader.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
ProofReader/ProofReader.o: Common/Sequence.h Common/PhaseSpace.h
ProofReader/ProofReader.o: /usr/include/stdio.h
ProofReader/ProofReader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
ProofReader/ProofReader.o: /usr/include/bits/types.h
ProofReader/ProofReader.o: /usr/include/bits/wordsize.h
ProofReader/ProofReader.o: /usr/include/bits/typesizes.h /usr/include/libio.h
ProofReader/ProofReader.o: /usr/include/_G_config.h /usr/include/wchar.h
ProofReader/ProofReader.o: /usr/include/bits/wchar.h /usr/include/gconv.h
ProofReader/ProofReader.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
ProofReader/ProofReader.o: /usr/include/bits/stdio_lim.h
ProofReader/ProofReader.o: /usr/include/bits/sys_errlist.h Common/PackedSeq.h
ProofReader/ProofReader.o: Common/HitRecord.h Common/SeqRecord.h
ProofReader/ProofReader.o: Common/CommonUtils.h
Partition/Partition.o: /usr/include/stdio.h /usr/include/features.h
Partition/Partition.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Partition/Partition.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
Partition/Partition.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
Partition/Partition.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Partition/Partition.o: /usr/include/_G_config.h /usr/include/wchar.h
Partition/Partition.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Partition/Partition.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
Partition/Partition.o: /usr/include/bits/stdio_lim.h
Partition/Partition.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
Partition/Partition.o: /usr/include/sys/stat.h /usr/include/bits/stat.h
Partition/Partition.o: /usr/include/math.h /usr/include/bits/huge_val.h
Partition/Partition.o: /usr/include/bits/mathdef.h
Partition/Partition.o: /usr/include/bits/mathcalls.h DataLayer/Reader.h
Partition/Partition.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Partition/Partition.o: /usr/include/assert.h Common/Sequence.h
Partition/Partition.o: Common/PhaseSpace.h Common/PackedSeq.h
Partition/Partition.o: Common/HitRecord.h Partition/Partition.h
Partition/Partition.o: DataLayer/IFileWriter.h DataLayer/PackedSeqWriter.h
Partition/Partition.o: Common/Config.h DataLayer/IFileReader.h
Partition/Partition.o: DataLayer/FastaReader.h DataLayer/PackedSeqReader.h
Partition/Partition.o: DataLayer/FastaWriter.h
Trimmer/Trimmer.o: Trimmer/Trimmer.h DataLayer/Reader.h Common/CommonDefs.h
Trimmer/Trimmer.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Trimmer/Trimmer.o: /usr/include/features.h /usr/include/sys/cdefs.h
Trimmer/Trimmer.o: /usr/include/gnu/stubs.h Common/Sequence.h
Trimmer/Trimmer.o: Common/PhaseSpace.h /usr/include/stdio.h
Trimmer/Trimmer.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
Trimmer/Trimmer.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
Trimmer/Trimmer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Trimmer/Trimmer.o: /usr/include/_G_config.h /usr/include/wchar.h
Trimmer/Trimmer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Trimmer/Trimmer.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
Trimmer/Trimmer.o: /usr/include/bits/stdio_lim.h
Trimmer/Trimmer.o: /usr/include/bits/sys_errlist.h Common/PackedSeq.h
Trimmer/Trimmer.o: Common/HitRecord.h Common/PartitionLoader.h
Trimmer/Trimmer.o: Common/Config.h DataLayer/PackedSeqWriter.h
Trimmer/Trimmer.o: DataLayer/IFileWriter.h
RemoveMultiplicity/RemoveMultiplicity.o: RemoveMultiplicity/RemoveMultiplicity.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/Config.h Common/CommonDefs.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/ReadPrb.h Common/Prb.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/assert.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/features.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/sys/cdefs.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/gnu/stubs.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/PhaseSpace.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/stdio.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stddef.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/types.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/wordsize.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/typesizes.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/libio.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/_G_config.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/wchar.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/wchar.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/gconv.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/lib/gcc/i386-redhat-linux/3.4.6/include/stdarg.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/stdio_lim.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/sys_errlist.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/Sequence.h Common/PackedSeq.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/HitRecord.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/Reader.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/PartitionLoader.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/PackedSeqWriter.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/IFileWriter.h
