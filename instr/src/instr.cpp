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

std::string MODE;

Instrument::Instrument(): ModulePass(ID), race(getenv("BENCH")) {
  MODE = getenv("MODE");
}

bool Instrument::runOnModule(Module &M) {
  this->AA = &getAnalysis<AliasAnalysis>();
  this->mod = &M;
  this->ctx = &getGlobalContext();

  for (Module::iterator F = M.begin(); F != M.end(); F ++) {
    if (F->getName() == "main") {
      // call init before main is about to start
      insert_before(F->begin()->begin(),
                    create_call(void_f("init"), NULL));
    }

    for (Function::iterator B = F->begin(); B != F->end(); B ++) {
      for (BasicBlock::iterator I = B->begin(); I != B->end(); I ++) {
        if (MODE == "NOTHING") continue;
        Instruction *inst = &*I;
        if (LoadInst *ld = dyn_cast<LoadInst>(inst)) {
          instr_load(ld);
        }
        if (StoreInst *st = dyn_cast<StoreInst>(inst)) {
          instr_store(st);
        }
        
        // not used in shared memory dependence tracing
        /*
        if (CallInst *call = dyn_cast<CallInst>(inst)) {
          instr_call(call);
        }
        */
      }
    }
  }

  auto *int_type = IntegerType::get(*ctx, 32);
  for (auto rm = remove_list.begin(); rm != remove_list.end(); rm ++) {
    (*rm)->replaceAllUsesWith(ConstantInt::get(int_type, 0));
    (*rm)->eraseFromParent();
  }
  remove_list.clear();

  return false; 
}

char Instrument::ID = 0;
static RegisterPass<Instrument> X("instrument", "Instrumentation", false, false);


