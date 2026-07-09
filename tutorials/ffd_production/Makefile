APP := build/bin/ffd_production
SRC := src/main.cc
OBJ := build/obj/main.o

UNAME_S := $(shell uname -s)
SYSTEM_OPENBLAS := $(firstword $(wildcard /usr/lib64/libopenblas.so.0 /lib64/libopenblas.so.0 /usr/lib/x86_64-linux-gnu/libopenblas.so.0))

CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra -pedantic

ifeq ($(UNAME_S),Darwin)
	CXX ?= clang++
	CXXFLAGS += -DACCELERATE_NEW_LAPACK
	LDFLAGS ?= -framework Accelerate
else
	ifneq ($(SYSTEM_OPENBLAS),)
		LDFLAGS ?= $(SYSTEM_OPENBLAS)
	else
		LDFLAGS ?= -llapack -lblas
	endif
endif

.PHONY: all clean print-config

all: $(APP)

$(APP): $(OBJ)
	@mkdir -p build/bin data
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(APP) $(LDFLAGS)

build/obj/main.o: $(SRC)
	@mkdir -p build/obj
	$(CXX) $(CXXFLAGS) -c $(SRC) -o build/obj/main.o

clean:
	rm -rf build

print-config:
	@echo "UNAME_S=$(UNAME_S)"
	@echo "CXX=$(CXX)"
	@echo "CXXFLAGS=$(CXXFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"
	@echo "SYSTEM_OPENBLAS=$(SYSTEM_OPENBLAS)"
