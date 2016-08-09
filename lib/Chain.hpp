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

#ifndef CHAIN_H
#define CHAIN_H

#include <vector>
#include <string>

//Forward Declaration
class Residue;
class Atom;

class Chain {
  private:
    std::vector<Residue *> resVec; //Coor of residue pointers
    std::vector<Atom *> atmVec; //Coor of atom pointers
//    bool sel;

  public:
    Chain();

    void reset();

    void addResidue(Residue* resEntry);
    void addAtom(Atom* atmEntry);

    //Get Chain info
    Atom* getAtom(const unsigned int& element);
    Residue* getResidue(const unsigned int& element);
    std::string getChainId();
    unsigned int getAtmVecSize();
    unsigned int getResVecSize();
//    void setSel(bool selin);
//    bool& getSel();
    void selAll();
    void deselAll();
};

#endif
