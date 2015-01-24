/*
 * RWTrace - program instrumentation
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

#include "instr.h"

void Instrument::instr_load(LoadInst *ld) {
  if (!is_potential_race(ld)) return;
  Value *ptr = ld->getPointerOperand();
  BitCastInst *cast;

  if (!need_read_retry()) {
    // need not to retry read
    ld->setVolatile(true);
    cast = create_cast(ptr);
    insert_before(ld, cast);
    insert_before(ld, create_call(void_f_ptr("before_read"), cast, NULL));
    insert_after(ld, create_call(void_f_ptr("after_read"), cast, NULL));
  } else {
    // need to retry read
    ld->setVolatile(false); // ld is useless
    unsigned alignment = ld->getAlignment();
    LoadInst *ld_clone = new LoadInst(ptr, "dynlk_ld_cl", true, alignment);
    LoadInst *ld_retry = new LoadInst(ptr, "dynlk_ld_rt", true, alignment);
    cast = create_cast(ptr);
    CallInst* read_judge = create_call(bool_f_ptr("after_tryread"), cast, NULL);
    SelectInst *verdict = SelectInst::Create(read_judge, ld_retry, ld_clone);

    Instruction *insts[] = {
      cast,
      create_call(void_f_ptr("before_tryread"), cast, NULL),
      ld_clone,                                  // make a clone of R(x)
      read_judge,                                // R(x) inter-procedural? 
      create_call(void_f_ptr("before_read"), cast, NULL), // call before_read
      ld_retry,                                  // read retry
      create_call(void_f_ptr("after_read"), cast, NULL), // call after_read
      verdict                                    // R(x) = read_judge ? ld_retry : ld_clone
    };

    for (int j = (sizeof(insts) / sizeof(Value*)) - 1; j >= 0; j --) {
      insert_after(ld, insts[j]);
    }

    ld->replaceAllUsesWith(verdict); // replace all uses of ld
  }
}

void Instrument::instr_store(StoreInst *st) {
  if (!is_potential_race(st)) return;
  st->setVolatile(true);
  BitCastInst *ptr = create_cast(st->getPointerOperand());
  insert_before(st, ptr);
  insert_before(st, create_call(void_f_ptr("before_write"), ptr, NULL));
  insert_after(st, create_call(void_f_ptr("after_write"), ptr, NULL));
}
