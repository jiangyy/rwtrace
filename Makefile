CXXFLAGS = $(shell llvm-config --cxxflags) -g -O2
LDFLAGS  = $(shell llvm-config --ldflags)

all: bin/runtime.bc bin/instr.so

INSTR_SRCS = $(shell find instr/src/ -name "*.cpp")
INSTR_HDRS = $(shell find instr/include/ -name "*.hpp")

%.so: $(INSTR_SRCS) $(INSTR_HDRS)
	clang++ -shared -std=c++11 -Iinstr/include/ $(CXXFLAGS) $(LDFLAGS) -lLLVM-3.4 $(INSTR_SRCS) -o $@

RUNTIME_SRCS = $(shell find runtime/src/ -name "*.cpp")
RUNTIME_OBJS = $(RUNTIME_SRCS:.cpp=.bc)

RT_HDRS = $(shell find runtime/include/ -name "*.h")

bin/runtime.bc: $(RUNTIME_OBJS)
	llvm-link -o $@ $(RUNTIME_OBJS)

%.bc: %.cpp $(RT_HDRS)
	clang -DDYNLK -c -std=c++11 -g -O2 -emit-llvm -Iruntime/include -o $@ $(@:.bc=.cpp)

clean:
	rm -f bin/*.so bin/*.bc $(RUNTIME_OBJS) a.out *.bc *.ll infile
