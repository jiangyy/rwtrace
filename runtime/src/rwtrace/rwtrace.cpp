/*
 * RWTrace: tracing shared memory dependences
 * Copyright (C) 2014  Yanyan Jiang <jiangyy@outlook.com>
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "runtime.h"
#include "rwtrace.h"
#include "log.h"
#include "hash.h"
#include <sys/syscall.h>
#include <sys/types.h>
#include <cassert>


static void _rwtrace_finalize();

// The lock of current shared memory access L(a)
static __thread Lock *clk;
// Does the read pass the thread-locality test?
static __thread bool need_retry;
// Lock owner cache's structure
static __thread LockOwnerCache tlc;
static __thread long last_reader;

static Lock*
get_lock(void *ptr) {
  return &locks[LK(ptr)];
}

extern "C" {

void
_RWTRACE_init() {
  atexit(_rwtrace_finalize);
  init_done = true;
//  thread_create_order.push_back(gettid());
  start_timer();
}

static void
_rwtrace_finalize() {
  fprintf(stderr, "[RWTrace] finished in %.2lfs.\n", get_time());
  flush_log(LOG_FILENAME);
  all_done();
}

void 
_RWTRACE_before_tryread(void *ptr) {
#ifdef TRACE_WAR
  if (!NEED_TRACK(ptr)) return;
  clk = get_lock(ptr);
  Hash *reader = clk->reader;
  if (reader != NULL) {
    last_reader = reader->quick_update(gettid(), evt_id);
  } else {
    last_reader = -1;
  }
  __sync_synchronize(); // MFENCE
#endif
}

bool
_RWTRACE_after_tryread(void *ptr) {
  if (!NEED_TRACK(ptr)) return false;
#ifndef TRACE_WAR
  clk = get_lock(ptr);
#endif
  // read versions AFTER read is performed
  //   there is no memory ordering issue because
  //   these fields are volatile
  bool  valid = clk->valid;
  pid_t tid = clk->last_writer;
  long  evid = clk->version;
  need_retry = ( !valid || !tlc.hit(clk, tid, evid) );
  return need_retry;
}

void
_RWTRACE_before_read(void *ptr) {
  if (!NEED_TRACK(ptr)) return;
  thread_local_event();

  if (need_retry) {
    clk = get_lock(ptr);
    clk->write_lock();
#ifdef TRACE_WAR
    Hash *hash = clk->reader;
    if (hash == NULL) {
      hash = hash_alloc();
      clk->reader = hash;
    } 
    hash->sync_update(gettid(), evt_id);
#endif
  }
}

void
_RWTRACE_after_read(void *ptr) {
  if (!NEED_TRACK(ptr)) return;
  if (need_retry) {
    pid_t clk_tid = clk->last_writer;
    long clk_ver = clk->version;
    clk->write_unlock();
    tlc.put(clk, clk_tid, clk_ver);
    LOG_RAW(clk_tid, clk_ver, evt_id);
    // There may be a over-written WAR dependence
    //   and [last_reader is the correct reader]
#ifdef TRACE_WAR
    if (last_reader != -1) {
      LOG_WAR(-1, evt_id, last_reader);
    }
#endif
  }
}

static __thread pid_t last_tid;
static __thread long  last_ver;
void
_RWTRACE_before_write(void *ptr) {
  if (!NEED_TRACK(ptr)) return;
  thread_local_event();

  clk = get_lock(ptr);
  clk->write_lock();

  // update version BEFORE write
  last_tid = clk->last_writer;
  last_ver = clk->version;
  clk->last_writer = gettid();
  clk->version = evt_id;
}

static __thread Hash *backup;

void
_RWTRACE_after_write(void *ptr) {
  if (!NEED_TRACK(ptr)) return;

#ifdef TRACE_WAR
  __sync_synchronize(); // MFENCE
#endif

  // whether this write thread-local?
  bool tl_write = tlc.hit(clk, last_tid, last_ver);

#ifdef TRACE_WAR
  pid_t ctid = gettid();
  Hash *old_hash = clk->reader, *new_hash;

  // backup is like a "double buffer"
  if (backup != NULL) {
    new_hash = backup;
    backup = NULL;
  } else {
    new_hash = hash_alloc();
  }
  new_hash->sync_update(ctid, -1);
  clk->reader = new_hash;
#endif
  clk->write_unlock();

  if (!tl_write) {
    LOG_WAW(last_tid, last_ver, evt_id);
  }
  tlc.put(clk, gettid(), evt_id);

#ifdef TRACE_WAR
  if (old_hash) {
    for (int i = 0; i < HASH_SIZE; i ++) {
      for (HashNode *nd = old_hash->data[i]; nd; nd = nd->next) {
        if (nd->key != ctid && nd->value != -1) {
          // nd->key == ctid: thread_local event
          // nd->value == -1: WAW dependence (we insert into hash table after write)
          LOG_WAR(nd->key, nd->value, evt_id);
        }
      }
    }
  }
  if (backup == NULL && old_hash != NULL) {
    backup = old_hash;
    old_hash->clear();
  }
#endif

}

}

