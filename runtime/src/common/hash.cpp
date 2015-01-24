/*
 * RWTrace wait-free hash table implementation
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

#include "hash.h"
#include "runtime.h"
#include <atomic>
#include <cstring>
#include <cassert>

Hash *
hash_alloc() {
  // THIS IS MEMORY LEAK.
  //   However, this is not much a problem because allocated Hash structures
  //   are recycled by the program (backup).
  //   We tried pooled memory allocation, it essentially has the same effect
  return new Hash();
}

// Each thread has a local pool of hash table nodes
// This pool is an approximate "read-copy-update"
static __thread HashNode *node_pool[NODE_POOL_SIZE];
static __thread int pool_f = 0, pool_r = 0;

static inline int
pool_next(int idx) {
  return (idx + 1) % NODE_POOL_SIZE;
}

static HashNode *
node_alloc(pid_t key, long value, HashNode *next) {
  if (pool_next(pool_f) == pool_r) {
    HashNode *ret = node_pool[pool_f];
    ret->key = key;      
    asm volatile("":::"memory");
    // Memory ordering of W(key) -> W(value) is important
    //   but MFENCE is not needed here.
    ret->value = value;
    ret->next = next;
    pool_f = pool_next(pool_f);
    __sync_synchronize();
    return ret;
  } else {
    return new HashNode(key, value, next);
  }
}

static void
node_free(HashNode *nd) {
  if ( pool_next(pool_f) == pool_r ) { // the ring buffer is full
    // TODO: leave it to a generational GC (RCU-like)
    //   this is also memory leak (but rare).
  } else {
    node_pool[pool_r] = nd;
    pool_r = pool_next(pool_r);
  }
}

long
Hash::quick_update(pid_t key, long value) {
  // The hash table lookup contains intended data race
  int h = ((int)key & (HASH_SIZE - 1));
  for (HashNode *nd = data[h]; nd; nd = nd->next) {
    if (nd->key == key) {
      long ret = nd->value;
      nd->value = value;
      __sync_synchronize(); // MFENCE
      // After thread t recycles a HashNode object, this object sticks
      //   to thread t, and if (nd->key != key), there is unexpected
      //   interleaving
      // Passing the following assertion suggests that our implementation is
      //   CORRECT. Check objdump result to see memory accesses and MFENCE
      //   instructions are not optimized out.
      assert(nd->key == key);

      return ret;
    }
  }
  return -1;
}

void
Hash::sync_update(pid_t key, long value) {
  // sync_update of M(a) is called only when L(a) is held
  int h = ((int)key & (HASH_SIZE - 1));
  HashNode *st = data[h];
  for (HashNode *nd = st; nd; nd = nd->next) {
    if (nd->key == key) {
      nd->value = value;
      return;
    }
  }
  data[h] = node_alloc(key, value, st);
}

void
Hash::clear() {
  // Traverse the hash table to reclaim memory.
  for (int i = 0; i < HASH_SIZE; i ++) {
    for (HashNode *nd = data[i]; nd; nd = nd->next) {
      node_free(nd);
    }
    data[i] = NULL;
  }
}
