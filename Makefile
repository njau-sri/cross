# make
# make win32
# make win64
# make macos

TARGET := cross

CXXFLAGS := -O2 -std=c++11 -Wall -DNDEBUG
LDFLAGS  := -s -static

OBJ := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

EXE := $(TARGET)

ifeq ($(MAKECMDGOALS),win32)
CXX := i686-w64-mingw32-g++
EXE := $(addsuffix .exe,$(TARGET))
endif

ifeq ($(MAKECMDGOALS),win64)
CXX := x86_64-w64-mingw32-g++
EXE := $(addsuffix .exe,$(TARGET))
endif

ifeq ($(MAKECMDGOALS),macos)
CXX := clang++
LDFLAGS := $(filter-out -s -static,$(LDFLAGS))
endif

.PHONY: all win32 win64 macos clean

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $?

win32: all

win64: all

macos: all

clean:
	$(RM) $(OBJ) $(TARGET) $(TARGET).exe
