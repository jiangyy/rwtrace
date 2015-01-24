/*
 * RWTrace logging header
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

#ifndef __LOG_H__
#define __LOG_H__

#include <functional>

#ifdef __APPLE__
typedef int pid_t;
#endif

// Log entry
struct LogEntry {
  pid_t thread;
  long event_id;
  long local_order;

  LogEntry() {}
  LogEntry(pid_t t, long a, long b) {
    thread = t;
    event_id = a;
    local_order = b;
  }
};

const int LogBlockSize = 1024;

// Per-thread chained append-only log
struct LogBlock {
  int size;
  LogEntry log[LogBlockSize];
  LogBlock *next;

  LogBlock() {
    size = 0;
    next = NULL;
  }

  void push(pid_t a, long b, long c) {
    log[size].thread = a;
    log[size].event_id = b;
    log[size].local_order = c;
    size ++;
  }
};

// Write (human-readable) log to a file
void flush_log(std::string fname);

#endif
