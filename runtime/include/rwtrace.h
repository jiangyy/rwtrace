/*
 * RWTrace
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

#ifndef __RWTRACE_H__
#define __RWTRACE_H__

#include "log.h"
#include "config.h"
#include "runtime.h"
#include <atomic>
#include <unordered_map>
#include <vector>
#include <string>

extern __thread LogBlock *current_log_block;
extern std::unordered_map<pid_t, LogBlock*> all_logs;
// extern std::vector<pid_t> thread_create_order;

// Logging dependences
void LOG_DEP(pid_t, long evt_id, long local_order);
#define LOG_RAW LOG_DEP
#define LOG_WAR LOG_DEP
#define LOG_WAW LOG_DEP

// A thread-local shadow version number slot
struct Slot {
  bool valid;
  pid_t pid;
  long ver;

  Slot();
  Slot(pid_t tid, long ver);
};

// A thread-local shadow version number
//   implemented as a direct cache
struct LockOwnerCache {
  Slot slots[NR_SLOT];
  LockOwnerCache();
  Slot *get_slot(Lock *lk);
  bool hit(Lock *lk, pid_t pid, long ver);  
  void put(Lock *lk, pid_t pid, long ver);
};

#endif
