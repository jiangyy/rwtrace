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
#include <sstream>
#include <cstdarg>

bool Instrument::need_read_retry() {
  if (MODE == "RWTRACE") {
    // special treatment for RWTrace, because it reads a shared address twice
    return true;
  }
  return false;
}

std::string Instrument::get_mode() {
  return MODE;
}

Constant *Instrument::void_f_ptr(std::string name) {
  std::string _name = "_" + get_mode() + "_" + std::string(name);
  return mod->getOrInsertFunction(_name, Type::getVoidTy(*ctx), Type::getInt8PtrTy(*ctx), NULL);
}

Constant *Instrument::bool_f_ptr(std::string name) {
  std::string _name = "_" + get_mode() + "_" + std::string(name);
  return mod->getOrInsertFunction(_name, Type::getInt1Ty(*ctx), Type::getInt8PtrTy(*ctx), NULL);
}

Constant *Instrument::void_f(std::string name) {
  std::string _name = "_" + get_mode() + "_" + std::string(name);
  return mod->getOrInsertFunction(_name, Type::getVoidTy(*ctx), NULL);
}

CallInst *Instrument::create_call(Constant *func, ...) {
  std::vector<Value*> args;
  va_list vlist;
  va_start(vlist, func);
  while (true) {
    Value *val = va_arg(vlist, Value*);
    if (val == NULL) break;
    args.push_back(val);
  }
  return CallInst::Create(func, args);
}

BitCastInst *Instrument::create_cast(Value *ptr) {
  return new BitCastInst(ptr, Type::getInt8PtrTy(*ctx));
}

std::string Instrument::get_debug_info(Instruction *i) {
  if (MDNode *N = i->getMetadata("dbg")) {
    DILocation Loc(N);
    std::stringstream ss;
    ss << Loc.getFilename().str() << ":" << Loc.getLineNumber();
    return ss.str();
  }
  return "None";
}

bool Instrument::name_check(Value *val) {
  return (strncmp(val->getName().str().c_str(), "dynlk", 5) == 0);
}

void Instrument::insert_before(Instruction *I, Instruction *i) {
  I->getParent()->getInstList().insert(I, i);
}

void Instrument::insert_after(Instruction *I, Instruction *i) {
  I->getParent()->getInstList().insertAfter(I, i);
}

void Instrument::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.addRequiredTransitive<AliasAnalysis>();
  AU.addPreserved<AliasAnalysis>();
}

bool Instrument::is_potential_race(Instruction *I) {
  if (name_check(I)) {
    return false;
  }

  // use race analysis to determine whether we need to instrument
  std::string debug = get_debug_info(I);
  if (debug != "None" && !race.is_race(debug)) {
    return false;
  }
  return true;
}
