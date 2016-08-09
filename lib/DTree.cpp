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

#include "DTree.hpp"

#include "Misc.hpp"

DTree::DTree(){
  root=NULL;
  classMap.clear();
}

DTree::~DTree(){
  delDTree();
}

void DTree::delDTree(DTreeNode *node){
  if (node != NULL){
    delDTree(node->left);
    delDTree(node->right);
    delete node;
  }
}

void DTree::addDTree(double key, DTreeNode *node, unsigned int index, std::string classin){
  if (node->left != NULL){
    addDTree(key, node->left, index, classin);
  }
  else{
    node->left=new DTreeNode;
    node->left->key_value=key;
    node->left->inx=index;
    node->left->left=NULL; //Sets left child of child node to NULL
    node->left->right=NULL; //Sets right child of child node to NULL
    node->left->leaf=classin; //The class of this node node
  }
  if (node->right != NULL){
    addDTree(key, node->right, index, classin);
  }
  else{
    node->right=new DTreeNode;
    node->right->key_value=key;
    node->right->inx=index;
    node->right->left=NULL; //Sets left child of child node to NULL
    node->right->right=NULL; //Sets right child of child node to NULL
    node->right->leaf=classin; //The class of this node node
  }
}

void DTree::genDTree(DTreeNode *&node, std::vector<std::string> &t, unsigned int &inx, std::string delim){
  //Pre-order string expected
  //Example: 0.0:0 1.0:1 2.0:2 3.0:3 A B 4.0:4 C D 5.0:5 E F 6.0:6 7.0:7 8.0:8 G H 9.0:9 I J 10.0:10 K L 
  //Value:Index or Class

  std::vector<std::string> s;
  unsigned int i;
  double d;

  if (delim.compare(0,1," " ) == 0){
    std::cerr << "Error: Delimiter cannot be whitespace!" << std::endl;
    return;
  }

  if (inx < t.size()){
    node=new DTreeNode;
    node->left=NULL;
    node->right=NULL;

    Misc::splitStr(Misc::trim(t.at(inx)), delim, s, false);
    if (s.size() == 1){
      //std::cerr << t.at(inx) << std::endl;
      node->leaf=s.at(0);
      return;
    }
    else if (s.size() == 2){
      //std::cerr << t.at(inx) << std::endl;
      std::stringstream(s.at(0)) >> d;
      node->key_value=d;
      std::stringstream(s.at(1)) >> i;
      node->inx=i; //Base zero
      this->genDTree(node->left, t, ++inx, delim);
      this->genDTree(node->right, t, ++inx, delim);
    }
    else{
      //I shouldn't be here
    }
  }
  else{
    //Do nothing
  }
}

std::string DTree::getDTreeClass(DTreeNode *node, const std::vector<double> &fin){
  if (node->left == NULL && node->right == NULL){
    //Return class
    if (this->getClassSize() != 0){
      return this->getClass(node->leaf);
    }
    else{
      return node->leaf;
    }
  }
  else if (node->left == NULL || node->right == NULL){
    if (node->left != NULL){
      //Go deeper
      return getDTreeClass(node->left, fin);
    }
    else{
      //Go deeper
      return getDTreeClass(node->right, fin);
    }
  }
  else{
    if (node->inx >= fin.size()){
      std::cerr << "Error: Missing feature (inx = " << node->inx << " )" << std::endl;
      return "";
    }
    else if (fin.at(node->inx) <= node->key_value){
      //Go deeper
      return getDTreeClass(node->left, fin);
    }
    else{
      //Go deeper
      return getDTreeClass(node->right, fin);
    }
  }
}

void DTree::addDTree(double key, unsigned int index, std::string classin){
  //Public version of the addDTree function
  //Takes care of when the root node has not been initialized
  if (root != NULL){
    addDTree(key, root, index, classin);
  }
  else{
    root=new DTreeNode;
    root->key_value=key;
    root->inx=index;
    root->left=NULL;
    root->right=NULL;
    root->leaf=classin;
  }
}

void DTree::genDTree(std::vector<std::string> &t, std::string delim){
  //Public version of the genDTree function
  //Pre-order string expected
  //Example: 0.0:0 1.0:1 2.0:2 3.0:3 A B 4.0:4 C D 5.0:5 E F 6.0:6 7.0:7 8.0:8 G H 9.0:9 I J 10.0:10 K L
  //Value:Index or Class

  unsigned int inx;
  std::vector<std::string> s;
  unsigned int i; 
  double d;

  inx=0;

  if (delim.compare(0,1," " ) == 0){
    std::cerr << "Error: Delimiter cannot be whitespace!" << std::endl;
    return;
  }

  if (inx < t.size()){
    root=new DTreeNode;
    root->left=NULL;
    root->right=NULL;
  
    Misc::splitStr(Misc::trim(t.at(inx)), delim, s, false);
    if (s.size() == 1){
      //std::cerr << t.at(inx) << std::endl;
      root->leaf=s.at(0);
      return;
    }
    else if (s.size() == 2){
      //std::cerr << t.at(inx) << std::endl;
      std::stringstream(s.at(0)) >> d;
      root->key_value=d;
      std::stringstream(s.at(1)) >> i;
      root->inx=i; //Base zero
      this->genDTree(root->left, t, ++inx, delim);
      this->genDTree(root->right, t, ++inx, delim);
    }
    else{
      //I shouldn't be here
    }
  }
  else{
    //Do nothing
  }
}

void DTree::delDTree(){
  //Public version of the delDTree function
  //Starts delete recursion starting at root node
  delDTree(root);
}

DTreeNode* DTree::getDTreeRoot(){
  if (root != NULL){
    return root;
  }
  else{
    return NULL;
  }
}

std::string DTree::getDTreeClass(const std::vector<double> &fin){
  //Public version of the getDTreeClass function
  if (root != NULL){
    return getDTreeClass(root, fin);
  }
  else{
    return "?"; 
  }
}

DTreeNode* DTreeNode::getDTreeNodeLeft(){
  return this->left;
}

DTreeNode* DTreeNode::getDTreeNodeRight(){
  return this->right;
}

unsigned int DTree::getClassSize(){
  return classMap.size();
}

void DTree::addClass(std::string valin, std::string classin){
  classMap.insert(std::pair<std::string,std::string>(valin, classin));  
}

void DTree::delClass(std::string classin){
  
}

std::string DTree::getClass(std::string valin){
  if (classMap.find(valin) != classMap.end()){
    return classMap.at(valin);
  }
  else{
    std::cerr << "Error: Unrecognized class mapping \"" << valin << "\"" << std::endl;
    return "";
  }
}


