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

#include "Prmtop.hpp"

#include "Misc.hpp"

#include <fstream>
#include <algorithm>

Prmtop::Prmtop (){
  mass.clear();
  charge.clear();
  atmtype.clear();
  //Hardcode HIS atom types with doubly deprotonated side chain nitrogens
  atmtype.insert(std::make_pair(std::make_pair("HIS","N"), "NH1"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","HN"), "H"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","CA"), "CT1"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","HA"), "HB1"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","CB"), "CT2"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","HB1"), "HA2"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","HB2"), "HA2"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","ND1"), "NR2"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","HD1"), "H"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","CG"), "CPH1"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","CE1"), "CPH2"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","HE1"), "HR1"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","NE2"), "NR2"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","CD2"), "CPH1"));
  atmtype.insert(std::make_pair(std::make_pair("HIS","HD2"), "HR3"));
}

void Prmtop::readTopology(const std::string& topin){
  std::ifstream topFile;
  std::istream* topinp;
  std::string line;
  std::vector<std::string> s;
  double m; //mass
  double c; //charge
  std::map<std::string, double> massRef;
  std::string resname;
  std::string word;

  resname.clear();

  if (topin.length() > 0){
    topFile.open(topin.c_str(), std::ios::in);
    topinp=&topFile;

    if (topinp->good() == false){
      std::cerr << "Warning: Topology file \"" << topin;
      std::cerr << " could not be read!" << std::endl;
    }

    while (topinp->good() && !(topinp->eof())){
      getline(*topinp, line);
      word=Misc::toupper(Misc::trim(line).substr(0,4));
      if (word.compare(0,4,"MASS") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 4){
          std::stringstream(s.at(3)) >> m;
          massRef.insert(std::make_pair(s.at(2), m));
          atomTypes.push_back(s.at(2));
        }
      }
      else if (word.compare(0,4,"RESI") == 0 || Misc::trim(line).compare(0,4,"PRES") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 3){
          resname=s.at(1);
        }
      }
      else if (word.compare(0,4,"ATOM") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 3){
          atmtype.insert(std::make_pair(std::make_pair(resname, s.at(1)), s.at(2)));
          mass.insert(std::make_pair(std::make_pair(resname, s.at(1)), massRef[s.at(2)]));
          std::stringstream(s.at(3)) >> c;
          charge.insert(std::make_pair(std::make_pair(resname, s.at(1)), c));
        }
      }
      else{
        //Do nothing
      }
    }
  }

  //Sort, remove duplicates, and resize atom type vector
  std::sort(atomTypes.begin(), atomTypes.end());
  std::vector<std::string>::iterator it=std::unique(atomTypes.begin(), atomTypes.end());
  atomTypes.resize(distance(atomTypes.begin(),it));

}


void Prmtop::readParameter(const std::string& prmin){

}


std::string Prmtop::getAtmType(const std::string& resnamein, const std::string& atmnamein, bool verboseFlag){
  if (atmtype.find(std::make_pair(resnamein, atmnamein)) != atmtype.end()){
    return  atmtype[std::make_pair(resnamein, atmnamein)];
  }
  else{
    if (verboseFlag == true){
      std::cerr << "Warning: Could not a find an atom type for atom " << resnamein;
      std::cerr << " " << atmnamein << " and was set to \"DUMB\"" << std::endl;
    }
    return "DUMB";
  }
}


double Prmtop::getMass(const std::string& resnamein, const std::string& atmnamein, bool verboseFlag){
  if (mass.find(std::make_pair(resnamein, atmnamein)) != mass.end()){
    return  mass[std::make_pair(resnamein, atmnamein)];
  }
  else{
    if (verboseFlag == true){
      std::cerr << "Warning: Could not a find mass for atom " << resnamein;
      std::cerr << " " << atmnamein << " and was set to 1.0" << std::endl;
    }
    return 1.0;
  }
}

double Prmtop::getCharge(const std::string& resnamein, const std::string& atmnamein, bool verboseFlag){
  if (charge.find(std::make_pair(resnamein, atmnamein)) != charge.end()){
    return charge[std::make_pair(resnamein, atmnamein)];
  }
  else{
    if (verboseFlag == true){
      std::cerr << "Warning: Could not find a charge for atom " << resnamein;
      std::cerr << " " << atmnamein << " and was set to 0.0." << std::endl;
    }
    return 0.0;
  }
}

std::vector<std::string> Prmtop::getAtomTypes(){
  return atomTypes;
}
