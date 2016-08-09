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

#include "Chain.hpp"

#include "Residue.hpp"
#include "Atom.hpp"

Chain::Chain (){
  resVec.clear();
  atmVec.clear();
//  sel=true;
}

void Chain::reset(){
  resVec.clear();
  atmVec.clear();
//  sel=true;
}

void Chain::addResidue(Residue* resEntry){
  if (resEntry->getResId()){
    resVec.push_back(resEntry);
  }
}

void Chain::addAtom(Atom* atmEntry){
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

Atom* Chain::getAtom (const unsigned int& element){
  if (element >= atmVec.size()){
    return NULL;
  }
  else{
    return atmVec.at(element);
  }
}

Residue* Chain::getResidue (const unsigned int& element){
  if (element >= resVec.size()){
    return NULL;
  }
  else{
    return resVec.at(element);
  }
}

std::string Chain::getChainId(){
  return this->getAtom(0)->getChainId();
}

unsigned int Chain::getAtmVecSize(){
  return atmVec.size();
}

unsigned int Chain::getResVecSize(){
  return resVec.size();
}

//void Chain::setSel(bool selin){
//  sel=selin;
//}

//bool& Chain::getSel(){
//  return sel;
//}

void Chain::selAll(){
//  sel=true;
  unsigned int i;
//  for (i=0; i< this->getResVecSize(); i++){
//    this->getResidue(i)->setSel(true);
//  }
  for (i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(true);
  }
}

void Chain::deselAll(){
//  sel=false;
  unsigned int i;
//  for (i=0; i< this->getResVecSize(); i++){
//    this->getResidue(i)->setSel(false);
//  }
  for (i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(false);
  }  
}
