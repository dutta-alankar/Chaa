# Chaa build
#
#   make            optimized build with HDF5 output support
#   make HDF5=0     build without the HDF5 dependency
#   make debug      unoptimized build with bounds checking

CHPL      ?= chpl
HDF5      ?= 1
BIN        = bin/chaa
SRC        = $(wildcard src/*.chpl) src/chaa_h5.h
CHPLFLAGS ?= --fast

# locate libhdf5: honour HDF5_LIBDIR, else try h5cc/brew/Debian paths
ifeq ($(HDF5),1)
  HDF5_LIBDIR ?= $(shell \
    if [ -n "$$HDF5_DIR" ]; then echo $$HDF5_DIR/lib; \
    elif command -v brew >/dev/null 2>&1 && brew --prefix hdf5 >/dev/null 2>&1; then \
      echo $$(brew --prefix hdf5)/lib; \
    elif [ -d /usr/lib/x86_64-linux-gnu/hdf5/serial ]; then \
      echo /usr/lib/x86_64-linux-gnu/hdf5/serial; \
    elif [ -d /usr/lib/aarch64-linux-gnu/hdf5/serial ]; then \
      echo /usr/lib/aarch64-linux-gnu/hdf5/serial; \
    fi)
  HDF5_FLAGS = -shdf5Enabled=true -lhdf5
  ifneq ($(strip $(HDF5_LIBDIR)),)
    HDF5_FLAGS += -L$(HDF5_LIBDIR)
  endif
else
  HDF5_FLAGS =
endif

all: $(BIN)

$(BIN): $(SRC)
	@mkdir -p bin
	$(CHPL) $(CHPLFLAGS) $(HDF5_FLAGS) --main-module Chaa $(SRC) -o $(BIN)

debug: CHPLFLAGS = --checks -g
debug: $(BIN)

clean:
	rm -rf bin

.PHONY: all debug clean
