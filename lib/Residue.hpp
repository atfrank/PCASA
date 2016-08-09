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

#ifndef RESIDUE_H
#define RESIDUE_H

#include <vector>
#include <string>

//Forward Declaration
class Atom;

class Residue {
  private:
    std::vector<Atom*> atmVec;
//    bool sel;

  public:
    Residue();

    void reset();
    int getResId();
    std::string getResName();
    std::string getChainId();
    Atom* getStart();
    Atom* getEnd();
    std::string getSegId();
    void addAtom(Atom* atmEntry);
    std::vector<Atom*>& getAtmVec();
    Atom* getAtom (const unsigned int &element);
    unsigned int getAtmVecSize();
//    void setSel(bool selin);
//    bool& getSel();
    void selAll();
    void deselAll();

    static std::string aa321(const std::string &aa);
    static std::string aa123(const std::string &aa);
};

#endif
