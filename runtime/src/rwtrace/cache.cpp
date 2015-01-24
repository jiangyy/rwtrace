/*
 * RWTrace: per-thread shadow memory implementation
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

#include "rwtrace.h"

Slot::Slot() {
  valid = false;
}

Slot::Slot(pid_t p, long v) {
  valid = true;
  pid = p; ver = v;
}


LockOwnerCache::LockOwnerCache() {
}

Slot *
LockOwnerCache::get_slot(Lock *lk) {
  int h = (unsigned long)lk % NR_SLOT;
  return &slots[h];
}

bool
LockOwnerCache::hit(Lock *lk, pid_t pid, long ver) {
  // the direct mapping is efficient enough for ~10K slots
  Slot *slot = get_slot(lk);
  if (!slot->valid) return false;
  pid_t known_pid = slot->pid;
  long known_ver = slot->ver;
  return known_pid == pid && known_ver == ver;
}

void
LockOwnerCache::put(Lock *lk, pid_t pid, long ver) {
  Slot *slot = get_slot(lk);
  slot->valid = true;
  slot->pid = pid;
  slot->ver = ver;
}
