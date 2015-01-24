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

#ifndef __INSTR_H__
#define	__INSTR_H__

#include <set>
#include <string>
#include <llvm/Analysis/AliasAnalysis.h>
#include <llvm/IR/Instructions.h>
#include <llvm/Pass.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/Constants.h>
#include <llvm/DebugInfo.h>
using namespace llvm;

extern std::string MODE;


/*
  Static data race analysis results from RELAY
  By default disabled.
*/
class RaceReport {
private:
  std::set<std::string> racy;
  bool exists = false;
  bool enabled = false;

public:
  RaceReport(const char *);
  bool is_race(std::string);
};


/*
  The class for program instrumentation.
  Instrument class reads MODE environment variable, and generates procedure
    calls before and after an interested event
  This is an LLVM optimization pass
*/
class Instrument: public ModulePass {
public:
  static char ID;
  Instrument();

private:
  Module *mod;
  AliasAnalysis *AA;
  LLVMContext *ctx;
  RaceReport race;
  std::set<Instruction*> remove_list;
  
  std::string get_mode();
  bool need_read_retry();
  Constant *void_f_ptr(std::string name);
  Constant *bool_f_ptr(std::string name);
  Constant *void_f(std::string name);
  CallInst *create_call(Constant *func, ...);
  BitCastInst *create_cast(Value *ptr);
  std::string get_debug_info(Instruction *i);
  bool name_check(Value *val);
  void insert_before(Instruction *I, Instruction *i);
  void insert_after(Instruction *I, Instruction *i);
  bool is_potential_race(Instruction *I);
  
  void instr_load(LoadInst *ld);
  void instr_store(StoreInst *st);
  void instr_call(CallInst *call);
  
public:
  virtual void getAnalysisUsage(AnalysisUsage &AU) const;
  virtual bool runOnModule(Module &M);
};

#endif
