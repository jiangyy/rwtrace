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

static struct FunctionDesc {
  const char *name; 
  bool remove_at_replay;
  int pointer_pos; // we are interested in a pointer parameter
} pthread_lib[] = {
  {"_create", false, -1},
  {"_join", false, -1},
  {"_lock", true, 0},
  {"_unlock", true, 0},
  {"_wait", true, -1},
  {"_signal", true, -1},
  {"_broadcast", true, -1}
};

void Instrument::instr_call(CallInst* call) {
  Function *f = call->getCalledFunction();
  if (f) {
    std::string func = f->getName().str();
    if (func.find("pthread") != std::string::npos) {
      int nfound = 0;
      for (int i = 0; i < sizeof(pthread_lib) / sizeof(FunctionDesc); i ++) {
        FunctionDesc *desc = &pthread_lib[i];
        if (func.find(desc->name) != std::string::npos) {
          nfound ++;
          std::string before = std::string("before") + desc->name;
          std::string after  = std::string("after") + desc->name;
          if (desc->pointer_pos != -1) {
            Value *arg = call->getArgOperand(desc->pointer_pos);
            BitCastInst *casted = create_cast(arg);
            insert_before(call, casted);
            insert_before(call, create_call(void_f_ptr(before), casted, NULL));
            insert_after(call, create_call(void_f_ptr(after), casted, NULL));
          } else {
            insert_before(call, create_call(void_f(before), NULL));
            insert_after(call, create_call(void_f(after), NULL));
          }
          if (get_mode() == "REPLAY" && desc->remove_at_replay) {
            remove_list.insert(call);
          }

        }
      }
      assert(nfound <= 1);
      if (nfound != 1) {
        fprintf(stderr, "WARN: call to %s has %d matches!\n", func.c_str(), nfound);
      }
    }
  }
}
