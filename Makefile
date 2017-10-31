CC = g++ 
##CC =icpc   
TARGET = bMatching
CFLAGS =-fopenmp -O3 
LDFLAGS =  $(CFLAGS)
OBJECTS = \
	mtxReader.o \
	bSuitor.o \
	bSuitorD.o \
	$(TARGET).o

all: $(TARGET) message


$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS) 

.cpp.o: 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS)

message:
	echo "Executable: $(TARGET) has been created"

