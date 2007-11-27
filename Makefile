# Makefile for the LGAP Assembly suite


#########################################################################
# DEFINES
COMMON_DIR = Common
DATALAYER_DIR = DataLayer
LGAP_DIR = LGAP
PROOFREAD_DIR = ProofReader

ALL_DIRS = $(COMMON_DIR) $(DATALAYER_DIR) $(LGAP_DIR) $(PROOFREAD_DIR)
 
COMMON_OBJ = $(patsubst %.cpp,%.o,$(wildcard $(COMMON_DIR)/*.cpp))
DATALAYER_OBJ = $(patsubst %.cpp,%.o,$(wildcard $(DATALAYER_DIR)/*.cpp))
LGAP_OBJ = $(patsubst %.cpp,%.o,$(wildcard $(LGAP_DIR)/*.cpp))
PROOFREAD_OBJ = $(patsubst %.cpp,%.o,$(wildcard $(PROOFREAD_DIR)/*.cpp))

ALL_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(LGAP_OBJ) $(PROOFREAD_OBJ)
LGAP_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(LGAP_OBJ)
PROOFREAD_OBJECTS = $(COMMON_OBJ) $(DATALAYER_OBJ) $(PROOFREAD_OBJ)

BASE_PATH = ..
BIN_PATH = $(BASE_PATH)/bin

#########################################################################
# FLAGS
CC =		g++

LIB_DIR =	
OPTIMIZE = -O2
#DEBUG= -g -D DEBUG -Wall -Wno-sign-compare
#PROFILE = -pg
OPTIONS =

INCLUDES =	-I$(COMMON_DIR) -I$(DATALAYER_DIR)

CPPFLAGS = $(OPTIMIZE) $(OPTIONS) $(DEBUG) $(PROFILE) $(INCLUDES)
LDFLAGS = $(PROFILE)
#########################################################################
# Rules

%.o: %.cpp
		$(CC) -c $(CPPFLAGS) $< -o $@
		
all:
	make LGAP ProofReader
	
LGAP:		$(LGAP_OBJECTS)
	$(CC) $(LDFLAGS) -o $(BIN_PATH)/LGAP $(LGAP_OBJECTS)
	
ProofReader:		$(PROOFREAD_OBJECTS)
	$(CC) $(LDFLAGS) -o $(BIN_PATH)/ProofReader $(PROOFREAD_OBJECTS) 
depend:
	makedepend *.cpp

clean:
	$(RM) $(BIN_PATH)/PathWalker $(BIN_PATH)/ProofReader $(ALL_OBJECTS)

