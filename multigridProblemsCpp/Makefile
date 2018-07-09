CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp CSR_Vector.cpp CSR_Matrix.cpp MG_Level.cpp MG_Structure.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=PA3

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o

distclean: clean
	rm $(EXECUTABLE)
