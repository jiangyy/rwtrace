#include "rwtrace.h"
#include "runtime.h"

__thread LogBlock *current_log_block = NULL;
std::unordered_map<pid_t, LogBlock*> all_logs;
// Used in replay, and traced by systemtap, not currently available
// std::vector<pid_t> thread_create_order;
static Lock logblock_lock;

// std::atomic<long> total_order(1);

void LOG(pid_t dep_tid, long dep_evid, long lid) {
  // log a dependence
  if (current_log_block == NULL) {
    logblock_lock.write_lock();
    current_log_block = new LogBlock();
    all_logs[gettid()] = current_log_block;
    logblock_lock.write_unlock();
  }
  if (current_log_block->size >= LogBlockSize) {
    LogBlock *old = current_log_block;
    current_log_block = new LogBlock();
    old->next = current_log_block;
  }
  current_log_block->push(dep_tid, dep_evid, lid);
}

void
LOG_DEP(pid_t dep_tid, long dep_evid, long lid) {
  if (dep_tid == -1 && dep_evid == -1) {
    // A null dependence
  } else {
    LOG(dep_tid, dep_evid, lid);
  }
}

// Not used for tracing shared memory dependences
// void
// LOG_SYN(long lid) {
//   long tick = total_order ++;
//   LOG(-1, tick, lid);
// }

void
flush_log(std::string fname) {
  long tot = 0;
  FILE *fp = fopen(fname.c_str(), "w");

//  fprintf(fp, "%ld", thread_create_order.size() - 1);
//  for (int i = 1; i < thread_create_order.size(); i ++) {
//    fprintf(fp, " %d", thread_create_order[i]);
//  }
//  fprintf(fp, "\n");

  long tot_tl = 0, tot_gl = 0;
  for (auto it = all_logs.begin(); it != all_logs.end(); it ++) {
    fprintf(fp, "%d ", it->first);
    int sz = 0;
    for (LogBlock *p = it->second; p; p = p->next) sz += p->size;
    long tl_evt = 0;
    fprintf(fp, "%d", sz);
    for (LogBlock *p = it->second; p; p = p->next) {
      for (int i = 0; i < p->size; i ++) {
        fprintf(fp, " %d,%ld,%ld", p->log[i].thread, p->log[i].event_id, p->log[i].local_order);
        tl_evt = std::max(tl_evt, p->log[i].local_order);
      }
    }
    tot_tl += tl_evt;
    tot_gl += sz;
    fprintf(fp, "\n");
  }
  fclose(fp);
  fprintf(stderr, "[RWTrace] %ld dependences logged out of %ld events.\n", tot_gl, tot_tl);
}
