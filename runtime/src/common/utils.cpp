/*
 * RWTrace common utility functions
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
#include <sys/syscall.h>
#include <sys/time.h>
#include <unistd.h>
#include <cstdio>
#include <signal.h>

__thread int __tid = -1;
__thread long evt_id = 1;
bool init_done = false;
pid_t first_thread = -1;

pid_t gettid() {
  // System calls are very slow. Therefore the thread id is buffered in the TLS
  if (__tid == -1) {
  #ifdef __APPLE__
    __tid = pthread_mach_thread_np(pthread_self());
  #else
    __tid = syscall(SYS_gettid);
  #endif
  }
  return __tid;
}

long
thread_local_event() {
  // evt_id is the "current" event in this thread
  return evt_id ++;
}

void
in_spin_loop() {
  #ifndef __APPLE__
    pthread_yield();
  #else
    asm volatile("PAUSE");
  #endif
}

void
all_done() {
  exit(0);
}

static struct timeval tpstart, tpend;

void
start_timer() {
  gettimeofday(&tpstart, NULL);
}

double
get_time() {
  gettimeofday(&tpend, NULL);
  return (tpend.tv_sec - tpstart.tv_sec) + (tpend.tv_usec - tpstart.tv_usec) / 1000000.0;
}
