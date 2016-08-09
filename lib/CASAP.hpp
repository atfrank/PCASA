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

#ifndef CASAP_H
#define CASAP_H

#include <string>

class CASAP {
  private:
    //Don't put std::vector<std::string> t here!
    //Can't instantiate object when using static call!

  public:
    static std::string getTree (unsigned int elem);
    static unsigned int getNTree();

};

#endif
