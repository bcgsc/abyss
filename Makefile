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
