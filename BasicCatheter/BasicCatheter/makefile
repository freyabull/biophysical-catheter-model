CXX = g++

BINARY = Catheter
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:%.cpp=%.o)


RM = rm -rf

all: $(BINARY)

%.o: %.cpp
	$(CXX) $(CPPFLAGS)  -o $@ -c $<

$(BINARY): $(OBJECTS)
	$(CXX) $(LFLAGS)  -o $@ $^  $(LIBS)

.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(BINARY)

