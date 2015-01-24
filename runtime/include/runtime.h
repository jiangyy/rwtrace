/*
 * RWTrace runtime information
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

#ifndef __RUNTIME_H__
#define __RUNTIME_H__

#include "config.h"
#include <pthread.h>
#include <atomic>
#include <map>
#include <unordered_map>
#include <vector>
#include "hash.h"
#ifdef __APPLE__
#include <sys/_types/_pid_t.h>
#endif

extern bool init_done;

// The base class of any lock, contians
struct LockInfo {
  // fields are volatile to ensure memory ordering
  pid_t volatile last_writer;
  bool  volatile valid;
  long  volatile version;

#ifdef TRACE_WAR
  // reader hash table is only useful in tracing WAR dependences
  Hash *reader;
#endif
};

struct MutexLock: public LockInfo {
  pthread_mutex_t lk;

  MutexLock();
  void read_lock();
  void read_unlock();
  void write_lock();
  void write_unlock();
};

struct ReadWriteLock: public LockInfo {
  std::atomic<int> num;
  pthread_mutex_t lk;

  ReadWriteLock();
  void read_lock();
  void read_unlock();
  void write_lock();
  void write_unlock();
};

// The runtime will use Lock only
#ifdef USE_MUTEX_LOCK
typedef MutexLock Lock;
#endif
#ifdef USE_RW_LOCK
typedef ReadWriteLock Lock;
#endif

// The lock pool
extern Lock locks[];

// return current thread's identifier
pid_t gettid();

// thread-local event counter
long thread_local_event(void);
extern __thread long evt_id;

// user-level hlt implmentation
extern void in_spin_loop();

// timing statistc
void start_timer();
double get_time();

// called when everything is okay
void all_done();

#endif
