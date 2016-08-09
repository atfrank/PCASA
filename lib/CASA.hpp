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

#ifndef CASA_H
#define CASA_H

#include <string>
#include <vector>
#include <map>


//#include "Molecule.hpp"

/* Forward Declaration  (only valid for pointers and references) */
class Molecule;

class CASA {
    private:
         std::vector<std::string> aminoAcids;
         std::map<std::string,double> refSASA;
         std::map<std::string,double> sasaCoef;
         std::map<std::string,double> hist;
         std::map<std::string,double> targetSASA; 
    public:
        CASA (Molecule *mol=NULL, const std::string fparmfile="", const std::string freffile="", const std::string fsasa="");
        void initializeRefSASA();
        void initializeSASACoef();
        void initializeAminoAcids();
        void initializeHist();
        void setHist(const std::string &key, const double &val);
        void accumHist(const std::string &key, const double &val);
        double predictHist(const std::string &resname);
        void printHist();
        void printHistHeader();
        void clearHist();
        double getHist(const std::string &key);
        double getRefSASA(const std::string &key);
        double getSASACoef(const std::string &key);
        double getTargetSASA(const std::string &key);
        void loadSASAFile(const std::string fsasa);
        void loadParmFile(const std::string fparmfile);
        void loadRefFile(const std::string freffile);
				std::string printRefSASAVec(const std::string &key);
				std::string printRefSASAVecHeader();
        
};
#endif
