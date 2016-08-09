//Sean M. Law
//Aaron T. Frank
    
/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SELECT_H
#define SELECT_H

#include <vector>
#include <map>
#include <string>

class Molecule;
class Atom;

class Select {
  private:
    std::map<std::string, std::string> selKeysAtm;
    std::map<std::string, std::string> selKeysRes;

  public:
    static void makeSel(Molecule* mol, std::string selin, bool dieFlag=true, bool verbose=true);
    void parseSel(std::string selin);

    //Recursive Descent Parser (RDP)
    std::vector<Atom*> recursiveDescentParser (const std::string &str, const std::vector<Atom *> &ref, const std::string &group="");
    static std::string getSelValue(const std::string &key);
    void initKeys(Molecule *mol);
    bool heavy(const std::string &str, const std::vector<std::string> &heavyVec);
    bool atom(const std::string &str, std::string typein, const std::vector<std::string> &AtomVec);
};

#endif
