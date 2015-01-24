## RWTrace
A shared memory dependence tracing tool for dynamic concurrent program analysis.

> RWTrace only traces shared memory dependence code. We plan to merge in more code (e.g., deterministic replay and predictive trace analysis tools) to this repository.

### Compiling the source

The source should compile using `make` if dependences are satisfied. We have tested RWTrace on 64bit Linux and LLVM3.4.

### Configuring RWTrace

The header file `runtime/include/config.h` contains some settings. Specifically, one can indicate whether writer-after-read dependences are traced (by defining or undefining `TRACE_WAR`). The dependence log file is outputed to `/dev/null` by default, which can be changed to any path to the system, or `/dev/stdout` for a screen output. Other parameters can be tuned (e.g., total amount of locks, and per-thread shadow memory size, etc.).

### Running RWTrace
Use the `run` script to trace shared memory dependences. This script contains essential commands to instrument the LLVM bitcode, and finally generate a runnable binary. 

    ./run water
    
runs the `water` benchmark (benchmarks are in `example/`).

### The dependence log
Each line of the dependence log is of the following form.

    tid n t1,reid1,leid1 t2,reid2,leid2 ...
    
`tid` denotes this thread's unique identifier. `n` denotes the amount of dependences. Each tuple `(t,reid,leid)` indicates a shared memory dependence of `(t,reid)`->`(tid,leid)`. A `t=-1` event is used to restore a potentially overwritten write-after-read dependence.