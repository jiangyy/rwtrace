/*
 * RWTrace lock implementations
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
#include <cassert>
#include <mutex>

Lock locks[NLOCKS] __attribute__((aligned(64)));

Lock global_lock; // this is not used in RWTrace

MutexLock::MutexLock() {
  pthread_mutex_init(&lk, 0);
  valid = true;
  last_writer = -1;
  version = -1;
}

void
MutexLock::read_lock() {
  pthread_mutex_lock(&lk);
}

void
MutexLock::read_unlock() {
  pthread_mutex_unlock(&lk);
}

void
MutexLock::write_lock() {
  pthread_mutex_lock(&lk);
}

void
MutexLock::write_unlock() {
  pthread_mutex_unlock(&lk);
  __sync_synchronize(); // MFENCE, to ensure sequential consistency
}


// Counter num is a std::atomic<int>
const int MAX_THREAD = 10000;
const int MAX_BACKOFF = 64;

int min(int a, int b) { // apple's clang does not recognize std::min
    return (a < b) ? a : b;
}

ReadWriteLock::ReadWriteLock() {
  num = MAX_THREAD;
  pthread_mutex_init(&lk, 0);
  last_writer = -1;
  valid = true;
  version = -1;
}

void
ReadWriteLock::read_lock() {
  int volatile count = -- num;
  int backoff = 1;
  while (count < 0) { // spin until there is no writer
    num ++;
    for (int i = 0; i < backoff; i ++) {
      in_spin_loop();
    }
    backoff = min(backoff * 2, MAX_BACKOFF);
    count = -- num;
  }
}

void
ReadWriteLock::read_unlock() {
  num ++;
}

void
ReadWriteLock::write_lock() {
  int volatile count = (num -= MAX_THREAD);
  int backoff = 1;
  while (count % MAX_THREAD != 0) {
    count = num.load();  // spin until no reader
    for (int i = 0; i < backoff; i ++) {
      in_spin_loop();
    }
    backoff = min(backoff * 2, MAX_BACKOFF);
  }
  pthread_mutex_lock(&lk);
}

void
ReadWriteLock::write_unlock() {
  num += MAX_THREAD;
  pthread_mutex_unlock(&lk);
}

