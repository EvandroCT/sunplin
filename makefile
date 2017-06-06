#in order to build device code, nvidia compiler must be invoked
CC=nvcc

#headers location
INC=include/

#sources location
SRC=src/

#created objects location
BUILD=build/

#libraries location
LIB=lib/

#list of objects to be created. TODO: automate process
OBJS=sunplin.o dtree.o htree.o phylosim.o cudutils.o

#-dc directive allows device code to be linked
CFLAGS=-w -std=c++11 -dc -I $(INC)

#linker flags.
LDFLAGS=-arch=sm_35

all: sunplin
	
sunplin: $(addprefix $(BUILD), $(OBJS))
	$(CC) $(LDFLAGS) $^ -o sunplin

$(BUILD)sunplin.o: $(SRC)sunplin.cu $(INC)htree.h $(INC)cudutils.h $(INC)phylosim.h
	$(CC) $(CFLAGS) $(SRC)sunplin.cu -o $(BUILD)sunplin.o

$(BUILD)dtree.o: $(SRC)dtree.cu $(INC)soatree.h
	$(CC) $(CFLAGS) $(SRC)dtree.cu -o $(BUILD)dtree.o

$(BUILD)htree.o: $(SRC)htree.cu $(INC)htree.h $(INC)cudutils.h
	$(CC) $(CFLAGS) $(SRC)htree.cu -o $(BUILD)htree.o

$(BUILD)phylosim.o: $(SRC)phylosim.cu $(INC)phylosim.h $(INC)cudutils.h
	$(CC) $(CFLAGS) $(SRC)phylosim.cu -o $(BUILD)phylosim.o

$(BUILD)cudutils.o: $(SRC)cudutils.cu
	$(CC) $(CFLAGS) $(SRC)cudutils.cu -o $(BUILD)cudutils.o

clean: 
	rm $(BUILD)*.o sunplin
