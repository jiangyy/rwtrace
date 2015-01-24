/*
 * RWTrace parameter settings
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

#ifndef __CONFIG_H__
#define __CONFIG_H__

// The lock implementation
//   (use pthread_mutex_lock by default, its a combined spin-mutex implementation)
#define USE_MUTEX_LOCK
//#define USE_RW_LOCK

// Are write-after-read dependences traced?
//#define TRACE_WAR

// NLOCKS locks, and each lock covers 2^(LKSCOPE) bytes, 
#define LKSCOPE  6
#define NLOCKS   (1 << 20)
// Shadow memory size
#define NR_SLOT 9997

// Mapping a pointer to its corresponding lock ID
#define LK(a) \
  (( ((unsigned long)(a)) >> LKSCOPE ) & (NLOCKS - 1))

// Initialization code and stack objects are not traced
#define NEED_TRACK(a) \
  ( init_done && ( (unsigned long)(a) < (unsigned long)0xffffffffff ) )

// Filename of the human-readable log
#define LOG_FILENAME "/dev/null"

#endif
