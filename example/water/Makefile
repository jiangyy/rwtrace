a.out: water.c
	clang -g -emit-llvm -c water.c -o a.bc
	opt -O2 < a.bc > opt.bc
	clang opt.bc -lpthread -lm
