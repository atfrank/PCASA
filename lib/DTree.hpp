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

#ifndef DTREE_H
#define DTREE_H

#include <string>
#include <vector>
#include <map>

struct DTreeNode{
  double key_value;
  unsigned int inx; //1-D vector<double> index, no response column!
  DTreeNode *left;
  DTreeNode *right;
  std::string leaf;

  DTreeNode* getDTreeNodeLeft();
  DTreeNode* getDTreeNodeRight();

};


class DTree {
  private:
    void delDTree(DTreeNode *node);
    DTreeNode* root;
    std::map<std::string, std::string> classMap;
    
    void addDTree(double key, DTreeNode *node, unsigned int index, std::string classin="");
    std::string getDTreeClass(DTreeNode *node, const std::vector<double> &fin);
    void genDTree(DTreeNode *&node, std::vector<std::string> &t, unsigned int &inx, std::string delim=":");

  public:
    DTree();
    ~DTree();

    void delDTree();

    void addDTree(double key, unsigned int index, std::string classin="");

    void genDTree(std::vector<std::string> &t, std::string delim=":");

    DTreeNode* getDTreeRoot();
    std::string getDTreeClass(const std::vector<double> &fin);

    unsigned int getClassSize();
    void addClass(std::string valin, std::string classin);
    void delClass(std::string classin);
    std::string getClass(std::string valin);
};

#endif


