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


#ifndef PDB_H
#define PDB_H

#include <map>
#include <string>

class Molecule;
class Residue;
class Atom;

class PDB {
  private:
    std::map<std::string, int> chnMap;
    std::string format; //Output format

  public:
    PDB();
    static void writePDBFormat (Molecule* mol, std::ostringstream &out, bool selFlag=true);
    static Molecule* readPDB (const std::string ifile, const int model=0, const std::string format="", const bool hetFlag=true, const bool remFlag=false);
    Atom* processAtomLine (std::string line, Atom* lastAtom);
    static std::string formatCHARMMResName (Atom* atmEntry);
    static int formatCHARMMResId(Atom* atmEntry, Residue* lastRes, Residue* nextRes);
};

#endif
