APP := build/bin/ffd_manual
SRC := src/main.cc
OBJ := build/obj/main.o

CXX := clang++
CXXFLAGS := -std=c++17 -O2 -DNDEBUG -Wall -Wextra -pedantic -DACCELERATE_NEW_LAPACK
LDFLAGS := -framework Accelerate

.PHONY: all clean

all: $(APP)

$(APP): $(OBJ)
	@mkdir -p build/bin
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(APP) $(LDFLAGS)

build/obj/main.o: $(SRC)
	@mkdir -p build/obj
	$(CXX) $(CXXFLAGS) -c $(SRC) -o build/obj/main.o

clean:
	rm -rf build
