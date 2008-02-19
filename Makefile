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
CC = g++
#CC =		mpic++

INCLUDE_DIR = $(HOME)/include
LIB_DIR = $(HOME)/lib
OPTIMIZE = -O2
64BIT = -m64
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
Common/CommonUtils.o: Common/CommonDefs.h /usr/include/stdlib.h
Common/CommonUtils.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/CommonUtils.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/CommonUtils.o: Common/SeqExt.h
Common/Config.o: Common/Config.h
Common/HitRecord.o: Common/HitRecord.h Common/CommonDefs.h
Common/HitRecord.o: /usr/include/stdlib.h /usr/include/features.h
Common/HitRecord.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/HitRecord.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/HitRecord.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/HitRecord.o: Common/PackedSeq.h Common/CommonUtils.h Common/Sequence.h
Common/HitRecord.o: Common/SeqExt.h
Common/PackedSeq.o: /usr/include/math.h /usr/include/features.h
Common/PackedSeq.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/PackedSeq.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
Common/PackedSeq.o: /usr/include/bits/mathcalls.h Common/PackedSeq.h
Common/PackedSeq.o: Common/CommonUtils.h Common/Sequence.h
Common/PackedSeq.o: Common/CommonDefs.h /usr/include/stdlib.h
Common/PackedSeq.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/PackedSeq.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/PackedSeq.o: Common/SeqExt.h
Common/PairRecord.o: Common/PairRecord.h Common/PackedSeq.h
Common/PairRecord.o: Common/CommonUtils.h Common/Sequence.h
Common/PairRecord.o: Common/CommonDefs.h /usr/include/stdlib.h
Common/PairRecord.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/PairRecord.o: /usr/include/gnu/stubs.h
Common/PairRecord.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/PairRecord.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/PairRecord.o: Common/SeqExt.h
Common/PartitionLoader.o: Common/PartitionLoader.h Common/CommonDefs.h
Common/PartitionLoader.o: /usr/include/stdlib.h /usr/include/features.h
Common/PartitionLoader.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/PartitionLoader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/PartitionLoader.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/PartitionLoader.o: Common/PhaseSpace.h /usr/include/stdio.h
Common/PartitionLoader.o: /usr/include/bits/types.h
Common/PartitionLoader.o: /usr/include/bits/wordsize.h
Common/PartitionLoader.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Common/PartitionLoader.o: /usr/include/_G_config.h /usr/include/wchar.h
Common/PartitionLoader.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Common/PartitionLoader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
Common/PartitionLoader.o: /usr/include/bits/stdio_lim.h
Common/PartitionLoader.o: /usr/include/bits/sys_errlist.h Common/Sequence.h
Common/PartitionLoader.o: Common/PackedSeq.h Common/CommonUtils.h
Common/PartitionLoader.o: Common/SeqExt.h Common/HitRecord.h Common/Config.h
Common/PhaseSpace.o: Common/PhaseSpace.h /usr/include/stdio.h
Common/PhaseSpace.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/PhaseSpace.o: /usr/include/gnu/stubs.h
Common/PhaseSpace.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/PhaseSpace.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
Common/PhaseSpace.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Common/PhaseSpace.o: /usr/include/_G_config.h /usr/include/wchar.h
Common/PhaseSpace.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Common/PhaseSpace.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
Common/PhaseSpace.o: /usr/include/bits/stdio_lim.h
Common/PhaseSpace.o: /usr/include/bits/sys_errlist.h Common/Sequence.h
Common/PhaseSpace.o: Common/CommonDefs.h /usr/include/stdlib.h
Common/PhaseSpace.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/PhaseSpace.o: Common/PackedSeq.h Common/CommonUtils.h Common/SeqExt.h
Common/PhaseSpace.o: Common/HitRecord.h
Common/Prb.o: /usr/include/math.h /usr/include/features.h
Common/Prb.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/Prb.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
Common/Prb.o: /usr/include/bits/mathcalls.h Common/Prb.h
Common/Prb.o: /usr/include/assert.h
Common/ReadPrb.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/ReadPrb.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/ReadPrb.o: /usr/include/gnu/stubs.h
Common/SeqExt.o: Common/SeqExt.h Common/CommonDefs.h /usr/include/stdlib.h
Common/SeqExt.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/SeqExt.o: /usr/include/gnu/stubs.h
Common/SeqExt.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/SeqExt.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/SeqRecord.o: Common/SeqRecord.h Common/PackedSeq.h
Common/SeqRecord.o: Common/CommonUtils.h Common/Sequence.h
Common/SeqRecord.o: Common/CommonDefs.h /usr/include/stdlib.h
Common/SeqRecord.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/SeqRecord.o: /usr/include/gnu/stubs.h
Common/SeqRecord.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/SeqRecord.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/SeqRecord.o: Common/SeqExt.h
Common/Sequence.o: /usr/include/stdio.h /usr/include/features.h
Common/Sequence.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Common/Sequence.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/Sequence.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
Common/Sequence.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Common/Sequence.o: /usr/include/_G_config.h /usr/include/wchar.h
Common/Sequence.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Common/Sequence.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
Common/Sequence.o: /usr/include/bits/stdio_lim.h
Common/Sequence.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
Common/Sequence.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
Common/Sequence.o: /usr/include/bits/mathcalls.h Common/Sequence.h
Common/Sequence.o: Common/CommonDefs.h /usr/include/stdlib.h Common/ReadPrb.h
Common/Sequence.o: Common/Prb.h /usr/include/assert.h
Common/SequencePair.o: Common/SequencePair.h Common/PackedSeq.h
Common/SequencePair.o: Common/CommonUtils.h Common/Sequence.h
Common/SequencePair.o: Common/CommonDefs.h /usr/include/stdlib.h
Common/SequencePair.o: /usr/include/features.h /usr/include/sys/cdefs.h
Common/SequencePair.o: /usr/include/gnu/stubs.h
Common/SequencePair.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/SequencePair.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Common/SequencePair.o: Common/SeqExt.h
Common/SimpleSequenceSpace.o: Common/SimpleSequenceSpace.h
Common/SimpleSequenceSpace.o: /usr/include/stdio.h /usr/include/features.h
Common/SimpleSequenceSpace.o: /usr/include/sys/cdefs.h
Common/SimpleSequenceSpace.o: /usr/include/gnu/stubs.h
Common/SimpleSequenceSpace.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Common/SimpleSequenceSpace.o: /usr/include/bits/types.h
Common/SimpleSequenceSpace.o: /usr/include/bits/wordsize.h
Common/SimpleSequenceSpace.o: /usr/include/bits/typesizes.h
Common/SimpleSequenceSpace.o: /usr/include/libio.h /usr/include/_G_config.h
Common/SimpleSequenceSpace.o: /usr/include/wchar.h /usr/include/bits/wchar.h
Common/SimpleSequenceSpace.o: /usr/include/gconv.h
Common/SimpleSequenceSpace.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
Common/SimpleSequenceSpace.o: /usr/include/bits/stdio_lim.h
Common/SimpleSequenceSpace.o: /usr/include/bits/sys_errlist.h
Common/SimpleSequenceSpace.o: Common/Sequence.h Common/CommonDefs.h
Common/SimpleSequenceSpace.o: /usr/include/stdlib.h Common/ReadPrb.h
Common/SimpleSequenceSpace.o: Common/Prb.h /usr/include/assert.h
Common/SimpleSequenceSpace.o: Common/PackedSeq.h Common/CommonUtils.h
Common/SimpleSequenceSpace.o: Common/SeqExt.h Common/HitRecord.h
DataLayer/FastaReader.o: DataLayer/FastaReader.h /usr/include/stdio.h
DataLayer/FastaReader.o: /usr/include/features.h /usr/include/sys/cdefs.h
DataLayer/FastaReader.o: /usr/include/gnu/stubs.h
DataLayer/FastaReader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
DataLayer/FastaReader.o: /usr/include/bits/types.h
DataLayer/FastaReader.o: /usr/include/bits/wordsize.h
DataLayer/FastaReader.o: /usr/include/bits/typesizes.h /usr/include/libio.h
DataLayer/FastaReader.o: /usr/include/_G_config.h /usr/include/wchar.h
DataLayer/FastaReader.o: /usr/include/bits/wchar.h /usr/include/gconv.h
DataLayer/FastaReader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
DataLayer/FastaReader.o: /usr/include/bits/stdio_lim.h
DataLayer/FastaReader.o: /usr/include/bits/sys_errlist.h
DataLayer/FastaReader.o: DataLayer/IFileReader.h Common/CommonDefs.h
DataLayer/FastaReader.o: /usr/include/stdlib.h Common/ReadPrb.h Common/Prb.h
DataLayer/FastaReader.o: /usr/include/assert.h Common/Sequence.h
DataLayer/FastaReader.o: Common/PackedSeq.h Common/CommonUtils.h
DataLayer/FastaReader.o: Common/SeqExt.h
DataLayer/FastaWriter.o: DataLayer/FastaWriter.h /usr/include/stdio.h
DataLayer/FastaWriter.o: /usr/include/features.h /usr/include/sys/cdefs.h
DataLayer/FastaWriter.o: /usr/include/gnu/stubs.h
DataLayer/FastaWriter.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
DataLayer/FastaWriter.o: /usr/include/bits/types.h
DataLayer/FastaWriter.o: /usr/include/bits/wordsize.h
DataLayer/FastaWriter.o: /usr/include/bits/typesizes.h /usr/include/libio.h
DataLayer/FastaWriter.o: /usr/include/_G_config.h /usr/include/wchar.h
DataLayer/FastaWriter.o: /usr/include/bits/wchar.h /usr/include/gconv.h
DataLayer/FastaWriter.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
DataLayer/FastaWriter.o: /usr/include/bits/stdio_lim.h
DataLayer/FastaWriter.o: /usr/include/bits/sys_errlist.h Common/CommonDefs.h
DataLayer/FastaWriter.o: /usr/include/stdlib.h Common/ReadPrb.h Common/Prb.h
DataLayer/FastaWriter.o: /usr/include/assert.h Common/Sequence.h
DataLayer/FastaWriter.o: DataLayer/IFileWriter.h Common/PackedSeq.h
DataLayer/FastaWriter.o: Common/CommonUtils.h Common/SeqExt.h
DataLayer/PackedSeqReader.o: DataLayer/PackedSeqReader.h /usr/include/stdio.h
DataLayer/PackedSeqReader.o: /usr/include/features.h /usr/include/sys/cdefs.h
DataLayer/PackedSeqReader.o: /usr/include/gnu/stubs.h
DataLayer/PackedSeqReader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
DataLayer/PackedSeqReader.o: /usr/include/bits/types.h
DataLayer/PackedSeqReader.o: /usr/include/bits/wordsize.h
DataLayer/PackedSeqReader.o: /usr/include/bits/typesizes.h
DataLayer/PackedSeqReader.o: /usr/include/libio.h /usr/include/_G_config.h
DataLayer/PackedSeqReader.o: /usr/include/wchar.h /usr/include/bits/wchar.h
DataLayer/PackedSeqReader.o: /usr/include/gconv.h
DataLayer/PackedSeqReader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
DataLayer/PackedSeqReader.o: /usr/include/bits/stdio_lim.h
DataLayer/PackedSeqReader.o: /usr/include/bits/sys_errlist.h
DataLayer/PackedSeqReader.o: DataLayer/IFileReader.h Common/CommonDefs.h
DataLayer/PackedSeqReader.o: /usr/include/stdlib.h Common/ReadPrb.h
DataLayer/PackedSeqReader.o: Common/Prb.h /usr/include/assert.h
DataLayer/PackedSeqReader.o: Common/Sequence.h Common/PackedSeq.h
DataLayer/PackedSeqReader.o: Common/CommonUtils.h Common/SeqExt.h
DataLayer/PackedSeqWriter.o: DataLayer/PackedSeqWriter.h /usr/include/stdio.h
DataLayer/PackedSeqWriter.o: /usr/include/features.h /usr/include/sys/cdefs.h
DataLayer/PackedSeqWriter.o: /usr/include/gnu/stubs.h
DataLayer/PackedSeqWriter.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/types.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/wordsize.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/typesizes.h
DataLayer/PackedSeqWriter.o: /usr/include/libio.h /usr/include/_G_config.h
DataLayer/PackedSeqWriter.o: /usr/include/wchar.h /usr/include/bits/wchar.h
DataLayer/PackedSeqWriter.o: /usr/include/gconv.h
DataLayer/PackedSeqWriter.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/stdio_lim.h
DataLayer/PackedSeqWriter.o: /usr/include/bits/sys_errlist.h
DataLayer/PackedSeqWriter.o: Common/CommonDefs.h /usr/include/stdlib.h
DataLayer/PackedSeqWriter.o: Common/ReadPrb.h Common/Prb.h
DataLayer/PackedSeqWriter.o: /usr/include/assert.h Common/Sequence.h
DataLayer/PackedSeqWriter.o: DataLayer/IFileWriter.h Common/PackedSeq.h
DataLayer/PackedSeqWriter.o: Common/CommonUtils.h Common/SeqExt.h
DataLayer/Reader.o: /usr/include/stdio.h /usr/include/features.h
DataLayer/Reader.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
DataLayer/Reader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
DataLayer/Reader.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
DataLayer/Reader.o: /usr/include/bits/typesizes.h /usr/include/libio.h
DataLayer/Reader.o: /usr/include/_G_config.h /usr/include/wchar.h
DataLayer/Reader.o: /usr/include/bits/wchar.h /usr/include/gconv.h
DataLayer/Reader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
DataLayer/Reader.o: /usr/include/bits/stdio_lim.h
DataLayer/Reader.o: /usr/include/bits/sys_errlist.h DataLayer/Reader.h
DataLayer/Reader.o: Common/CommonDefs.h /usr/include/stdlib.h
DataLayer/Reader.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
DataLayer/Reader.o: Common/Sequence.h Common/PhaseSpace.h Common/PackedSeq.h
DataLayer/Reader.o: Common/CommonUtils.h Common/SeqExt.h Common/HitRecord.h
DataLayer/Writer.o: DataLayer/Writer.h
ABYSS/Abyss.o: /usr/include/stdio.h /usr/include/features.h
ABYSS/Abyss.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
ABYSS/Abyss.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
ABYSS/Abyss.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
ABYSS/Abyss.o: /usr/include/bits/typesizes.h /usr/include/libio.h
ABYSS/Abyss.o: /usr/include/_G_config.h /usr/include/wchar.h
ABYSS/Abyss.o: /usr/include/bits/wchar.h /usr/include/gconv.h
ABYSS/Abyss.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
ABYSS/Abyss.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
ABYSS/Abyss.o: ABYSS/Abyss.h Common/Sequence.h Common/CommonDefs.h
ABYSS/Abyss.o: /usr/include/stdlib.h Common/ReadPrb.h Common/Prb.h
ABYSS/Abyss.o: /usr/include/assert.h DataLayer/Reader.h Common/PhaseSpace.h
ABYSS/Abyss.o: Common/PackedSeq.h Common/CommonUtils.h Common/SeqExt.h
ABYSS/Abyss.o: Common/HitRecord.h ABYSS/PathDriver.h ABYSS/Path.h
ABYSS/Abyss.o: Common/SeqRecord.h Common/PairRecord.h DataLayer/Writer.h
ABYSS/Abyss.o: Common/SequencePair.h Common/Config.h Common/PartitionLoader.h
ABYSS/Abyss.o: DataLayer/FastaWriter.h DataLayer/IFileWriter.h
ABYSS/Abyss.o: DataLayer/FastaReader.h DataLayer/IFileReader.h
ABYSS/Abyss.o: DataLayer/PackedSeqWriter.h Common/SimpleSequenceSpace.h
ABYSS/Path.o: ABYSS/Path.h Common/CommonDefs.h /usr/include/stdlib.h
ABYSS/Path.o: /usr/include/features.h /usr/include/sys/cdefs.h
ABYSS/Path.o: /usr/include/gnu/stubs.h
ABYSS/Path.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
ABYSS/Path.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
ABYSS/Path.o: Common/SeqRecord.h Common/PackedSeq.h Common/CommonUtils.h
ABYSS/Path.o: Common/Sequence.h Common/SeqExt.h
ABYSS/PathDriver.o: ABYSS/PathDriver.h Common/CommonDefs.h
ABYSS/PathDriver.o: /usr/include/stdlib.h /usr/include/features.h
ABYSS/PathDriver.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
ABYSS/PathDriver.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
ABYSS/PathDriver.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
ABYSS/PathDriver.o: ABYSS/Path.h Common/SeqRecord.h Common/PackedSeq.h
ABYSS/PathDriver.o: Common/CommonUtils.h Common/Sequence.h Common/SeqExt.h
ABYSS/PathDriver.o: Common/HitRecord.h Common/PhaseSpace.h
ABYSS/PathDriver.o: /usr/include/stdio.h /usr/include/bits/types.h
ABYSS/PathDriver.o: /usr/include/bits/wordsize.h
ABYSS/PathDriver.o: /usr/include/bits/typesizes.h /usr/include/libio.h
ABYSS/PathDriver.o: /usr/include/_G_config.h /usr/include/wchar.h
ABYSS/PathDriver.o: /usr/include/bits/wchar.h /usr/include/gconv.h
ABYSS/PathDriver.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
ABYSS/PathDriver.o: /usr/include/bits/stdio_lim.h
ABYSS/PathDriver.o: /usr/include/bits/sys_errlist.h Common/PairRecord.h
ABYSS/PathDriver.o: DataLayer/Writer.h Common/SequencePair.h
ProofReader/ProofReader.o: ProofReader/ProofReader.h DataLayer/Reader.h
ProofReader/ProofReader.o: Common/CommonDefs.h /usr/include/stdlib.h
ProofReader/ProofReader.o: /usr/include/features.h /usr/include/sys/cdefs.h
ProofReader/ProofReader.o: /usr/include/gnu/stubs.h
ProofReader/ProofReader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
ProofReader/ProofReader.o: Common/ReadPrb.h Common/Prb.h
ProofReader/ProofReader.o: /usr/include/assert.h Common/Sequence.h
ProofReader/ProofReader.o: Common/PhaseSpace.h /usr/include/stdio.h
ProofReader/ProofReader.o: /usr/include/bits/types.h
ProofReader/ProofReader.o: /usr/include/bits/wordsize.h
ProofReader/ProofReader.o: /usr/include/bits/typesizes.h /usr/include/libio.h
ProofReader/ProofReader.o: /usr/include/_G_config.h /usr/include/wchar.h
ProofReader/ProofReader.o: /usr/include/bits/wchar.h /usr/include/gconv.h
ProofReader/ProofReader.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
ProofReader/ProofReader.o: /usr/include/bits/stdio_lim.h
ProofReader/ProofReader.o: /usr/include/bits/sys_errlist.h Common/PackedSeq.h
ProofReader/ProofReader.o: Common/CommonUtils.h Common/SeqExt.h
ProofReader/ProofReader.o: Common/HitRecord.h Common/SeqRecord.h
Partition/Partition.o: /usr/include/stdio.h /usr/include/features.h
Partition/Partition.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Partition/Partition.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Partition/Partition.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
Partition/Partition.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Partition/Partition.o: /usr/include/_G_config.h /usr/include/wchar.h
Partition/Partition.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Partition/Partition.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
Partition/Partition.o: /usr/include/bits/stdio_lim.h
Partition/Partition.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
Partition/Partition.o: /usr/include/sys/stat.h /usr/include/bits/stat.h
Partition/Partition.o: /usr/include/math.h /usr/include/bits/huge_val.h
Partition/Partition.o: /usr/include/bits/mathdef.h
Partition/Partition.o: /usr/include/bits/mathcalls.h DataLayer/Reader.h
Partition/Partition.o: Common/CommonDefs.h Common/ReadPrb.h Common/Prb.h
Partition/Partition.o: /usr/include/assert.h Common/Sequence.h
Partition/Partition.o: Common/PhaseSpace.h Common/PackedSeq.h
Partition/Partition.o: Common/CommonUtils.h Common/SeqExt.h
Partition/Partition.o: Common/HitRecord.h Partition/Partition.h
Partition/Partition.o: DataLayer/IFileWriter.h DataLayer/PackedSeqWriter.h
Partition/Partition.o: Common/Config.h DataLayer/IFileReader.h
Partition/Partition.o: DataLayer/FastaReader.h DataLayer/PackedSeqReader.h
Partition/Partition.o: DataLayer/FastaWriter.h
Trimmer/Trimmer.o: Trimmer/Trimmer.h DataLayer/Reader.h Common/CommonDefs.h
Trimmer/Trimmer.o: /usr/include/stdlib.h /usr/include/features.h
Trimmer/Trimmer.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Trimmer/Trimmer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Trimmer/Trimmer.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Trimmer/Trimmer.o: Common/Sequence.h Common/PhaseSpace.h /usr/include/stdio.h
Trimmer/Trimmer.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
Trimmer/Trimmer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Trimmer/Trimmer.o: /usr/include/_G_config.h /usr/include/wchar.h
Trimmer/Trimmer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Trimmer/Trimmer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
Trimmer/Trimmer.o: /usr/include/bits/stdio_lim.h
Trimmer/Trimmer.o: /usr/include/bits/sys_errlist.h Common/PackedSeq.h
Trimmer/Trimmer.o: Common/CommonUtils.h Common/SeqExt.h Common/HitRecord.h
Trimmer/Trimmer.o: Common/PartitionLoader.h Common/Config.h
Trimmer/Trimmer.o: DataLayer/PackedSeqWriter.h DataLayer/IFileWriter.h
RemoveMultiplicity/RemoveMultiplicity.o: RemoveMultiplicity/RemoveMultiplicity.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/Config.h Common/CommonDefs.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/stdlib.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/features.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/sys/cdefs.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/gnu/stubs.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/ReadPrb.h Common/Prb.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/assert.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/PhaseSpace.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/stdio.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/types.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/wordsize.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/typesizes.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/libio.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/_G_config.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/wchar.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/wchar.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/gconv.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/stdio_lim.h
RemoveMultiplicity/RemoveMultiplicity.o: /usr/include/bits/sys_errlist.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/Sequence.h Common/PackedSeq.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/CommonUtils.h Common/SeqExt.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/HitRecord.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/Reader.h
RemoveMultiplicity/RemoveMultiplicity.o: Common/PartitionLoader.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/PackedSeqWriter.h
RemoveMultiplicity/RemoveMultiplicity.o: DataLayer/IFileWriter.h
Parallel/CommLayer.o: Parallel/CommLayer.h Common/PackedSeq.h
Parallel/CommLayer.o: Common/CommonUtils.h Common/Sequence.h
Parallel/CommLayer.o: Common/CommonDefs.h /usr/include/stdlib.h
Parallel/CommLayer.o: /usr/include/features.h /usr/include/sys/cdefs.h
Parallel/CommLayer.o: /usr/include/gnu/stubs.h
Parallel/CommLayer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Parallel/CommLayer.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Parallel/CommLayer.o: Common/SeqExt.h
Parallel/DistributedPhaseSpace.o: Parallel/DistributedPhaseSpace.h
Parallel/DistributedPhaseSpace.o: Common/SimpleSequenceSpace.h
Parallel/DistributedPhaseSpace.o: /usr/include/stdio.h
Parallel/DistributedPhaseSpace.o: /usr/include/features.h
Parallel/DistributedPhaseSpace.o: /usr/include/sys/cdefs.h
Parallel/DistributedPhaseSpace.o: /usr/include/gnu/stubs.h
Parallel/DistributedPhaseSpace.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Parallel/DistributedPhaseSpace.o: /usr/include/bits/types.h
Parallel/DistributedPhaseSpace.o: /usr/include/bits/wordsize.h
Parallel/DistributedPhaseSpace.o: /usr/include/bits/typesizes.h
Parallel/DistributedPhaseSpace.o: /usr/include/libio.h
Parallel/DistributedPhaseSpace.o: /usr/include/_G_config.h
Parallel/DistributedPhaseSpace.o: /usr/include/wchar.h
Parallel/DistributedPhaseSpace.o: /usr/include/bits/wchar.h
Parallel/DistributedPhaseSpace.o: /usr/include/gconv.h
Parallel/DistributedPhaseSpace.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
Parallel/DistributedPhaseSpace.o: /usr/include/bits/stdio_lim.h
Parallel/DistributedPhaseSpace.o: /usr/include/bits/sys_errlist.h
Parallel/DistributedPhaseSpace.o: Common/Sequence.h Common/CommonDefs.h
Parallel/DistributedPhaseSpace.o: /usr/include/stdlib.h Common/ReadPrb.h
Parallel/DistributedPhaseSpace.o: Common/Prb.h /usr/include/assert.h
Parallel/DistributedPhaseSpace.o: Common/PackedSeq.h Common/CommonUtils.h
Parallel/DistributedPhaseSpace.o: Common/SeqExt.h Common/HitRecord.h
Parallel/DistributedPhaseSpace.o: Parallel/CommLayer.h
Parallel/parallelAbyss.o: /usr/include/stdio.h /usr/include/features.h
Parallel/parallelAbyss.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
Parallel/parallelAbyss.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
Parallel/parallelAbyss.o: /usr/include/bits/types.h
Parallel/parallelAbyss.o: /usr/include/bits/wordsize.h
Parallel/parallelAbyss.o: /usr/include/bits/typesizes.h /usr/include/libio.h
Parallel/parallelAbyss.o: /usr/include/_G_config.h /usr/include/wchar.h
Parallel/parallelAbyss.o: /usr/include/bits/wchar.h /usr/include/gconv.h
Parallel/parallelAbyss.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
Parallel/parallelAbyss.o: /usr/include/bits/stdio_lim.h
Parallel/parallelAbyss.o: /usr/include/bits/sys_errlist.h
Parallel/parallelAbyss.o: Parallel/parallelAbyss.h Common/PackedSeq.h
Parallel/parallelAbyss.o: Common/CommonUtils.h Common/Sequence.h
Parallel/parallelAbyss.o: Common/CommonDefs.h /usr/include/stdlib.h
Parallel/parallelAbyss.o: Common/ReadPrb.h Common/Prb.h /usr/include/assert.h
Parallel/parallelAbyss.o: Common/SeqExt.h DataLayer/FastaReader.h
Parallel/parallelAbyss.o: DataLayer/IFileReader.h
Parallel/parallelAbyss.o: Parallel/DistributedPhaseSpace.h
Parallel/parallelAbyss.o: Common/SimpleSequenceSpace.h Common/HitRecord.h
Parallel/parallelAbyss.o: Parallel/CommLayer.h
