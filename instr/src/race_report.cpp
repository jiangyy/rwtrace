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
#include <string>
#include <fstream>
#include <regex.h>

RaceReport::RaceReport(const char *benchname) {
  if (!enabled) {
    return;
  }
  static char filename[1024];
  sprintf(filename, "%s/ciltrees/warnings.xml", benchname);
  static char buf[4096];

  regex_t re;
  regmatch_t value[3];
  // parsing RELAY's race report
  regcomp(&re, "file=\"([^\"]+)\" line=\"([^\"]+)\"", REG_EXTENDED);

  std::ifstream fin(filename);
  if (!fin.is_open()) {
    fprintf(stderr, "Cannot open RELAY analysis result file.\n");
    exists = false;
    return;
  }
  exists = true;

  while (!fin.eof()) {
    fin.getline(buf, 4096);
    if (regexec(&re, buf, 3, value, 0) != REG_NOMATCH) {
      std::string m1(buf + value[1].rm_so, buf + value[1].rm_eo);
      std::string m2(buf + value[2].rm_so, buf + value[2].rm_eo);
      std::string res = m1 + ":" + m2;

      racy.insert(res); // find a pair of racy statements
    }  
  }

  fin.close();
}

bool RaceReport::is_race(std::string a) {
  if (!enabled) return true;
  if (!exists) return true;
  // RELAY is expected to be sound
  return racy.find(a) != racy.end();
}
