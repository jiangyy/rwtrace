/*
 * RWTrace hash table header
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

#ifndef __HASH_H__
#define __HASH_H__

#include <unordered_map>

const int HASH_SIZE = 16;
const int NODE_POOL_SIZE = (1 << 10);
const int HASH_POOL_SIZE = (1 << 16);

// A <key, value> pair
struct HashNode {
  pid_t key;
  long value;
  HashNode *next;

  HashNode() {
  }

  HashNode(pid_t k, long v, HashNode *n) {
    key = k;
    value = v;
    next = n;
  }
};

// A set of <key, value> pairs
struct Hash {
  HashNode *data[HASH_SIZE];

  // update <key, value> pair, do nothing if lookup fails
  long quick_update(pid_t key, long value);
  // update <key, value> pair, insert it if lookup fails
  void sync_update(pid_t key, long value);
  // reclaim memory
  void clear();
};

Hash *hash_alloc();

#endif
