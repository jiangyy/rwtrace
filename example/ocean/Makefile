a.out: ocean.c
	clang -g -emit-llvm -c ocean.c -o a.bc
	opt -O2 < a.bc > opt.bc
	clang opt.bc -lpthread -lm
