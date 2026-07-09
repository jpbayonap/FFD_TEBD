PROJECT_ROOT := $(abspath $(CURDIR)/../..)
ITENSOR_DIR ?= $(PROJECT_ROOT)/vendor/ITensor

include $(ITENSOR_DIR)/this_dir.mk
include $(ITENSOR_DIR)/options.mk

TENSOR_HEADERS := $(ITENSOR_DIR)/itensor/all.h
SRC_DIR := $(CURDIR)/src
OBJ_DIR := $(CURDIR)/build/obj
BIN_DIR := $(CURDIR)/build/bin

OBJS := $(patsubst $(SRC_DIR)/%.cc,$(OBJ_DIR)/%.o,$(SRC))

.PHONY: all clean

all: $(BIN_DIR)/$(APP)

$(BIN_DIR)/$(APP): $(OBJS) $(ITENSOR_LIBS) $(TENSOR_HEADERS)
	@mkdir -p "$(BIN_DIR)"
	$(CCCOM) $(CCFLAGS) $(OBJS) -o "$@" $(LIBFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(TENSOR_HEADERS)
	@mkdir -p "$(OBJ_DIR)"
	$(CCCOM) -c $(CCFLAGS) "$<" -o "$@"

clean:
	rm -rf "$(CURDIR)/build"
