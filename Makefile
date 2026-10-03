# dnaln — build with `make` (release), `make debug` (ASan/UBSan), `make test`, `make bench`, `make demo`.
# Install: `make install` (PREFIX=/usr/local) or `make install PREFIX="$HOME/.local"`.
# Redistributable binary: `make PORTABLE=1` (baseline CPU features instead of -mcpu=native).
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(PORTABLE),1)
  # Binaries for redistribution: baseline CPU features only (SIMD still on: SSE4.1 / NEON).
  ifeq ($(UNAME_M),x86_64)
    ARCH ?= -march=x86-64-v2
  else ifeq ($(UNAME_S),Darwin)
    ARCH ?= -mcpu=apple-m1
  else
    ARCH ?= -march=armv8-a
  endif
  ifeq ($(UNAME_S),Linux)
    # Fully static (libc, libstdc++, zlib), so the binary has no runtime dependency on the
    # system's glibc or other shared libraries.
    LDFLAGS += -static
    LDLIBS  := -lz
  endif
  ifeq ($(UNAME_S),Darwin)
    # Runs on macOS 11 (Big Sur, the first release for Apple Silicon) and later.
    CXXFLAGS += -mmacosx-version-min=11.0
  endif
else ifeq ($(UNAME_M),x86_64)
  ARCH ?= -march=native
else
  ARCH ?= -mcpu=native
endif
OPT      ?= -O3 -DNDEBUG
CXXFLAGS += -std=c++17 $(OPT) $(ARCH) -Wall -Wextra -Wno-unused-parameter -MMD -MP
LDLIBS   ?= -lz
PREFIX   ?= /usr/local
BUILD    ?= build

LIB_SRCS := src/common.cpp src/seqio.cpp src/reference.cpp src/minimizer.cpp src/index.cpp \
            src/chain.cpp src/align.cpp src/align_simd.cpp src/mapper.cpp src/simulate.cpp \
            src/evaluate.cpp src/groundtruth.cpp src/bench.cpp src/demo.cpp
LIB_OBJS := $(LIB_SRCS:src/%.cpp=$(BUILD)/obj/%.o)

all: $(BUILD)/dnaln $(BUILD)/dnaln_tests

$(BUILD)/dnaln: $(LIB_OBJS) $(BUILD)/obj/main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/dnaln_tests: $(LIB_OBJS) $(BUILD)/obj/test_main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD)/obj/%.o: src/%.cpp | $(BUILD)/obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/obj/test_main.o: tests/test_main.cpp | $(BUILD)/obj
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

$(BUILD)/obj:
	mkdir -p $@

test: $(BUILD)/dnaln_tests
	./$(BUILD)/dnaln_tests

debug:
	$(MAKE) BUILD=build-debug OPT="-O0 -g -fsanitize=address,undefined -fno-omit-frame-pointer"

bench: $(BUILD)/dnaln
	./scripts/run_benchmark.sh

demo: $(BUILD)/dnaln
	./demo.sh

install: $(BUILD)/dnaln
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 $(BUILD)/dnaln "$(DESTDIR)$(PREFIX)/bin/dnaln"

uninstall:
	rm -f "$(DESTDIR)$(PREFIX)/bin/dnaln"

clean:
	rm -rf build build-debug

-include $(wildcard $(BUILD)/obj/*.d)
.PHONY: all test debug bench demo install uninstall clean
