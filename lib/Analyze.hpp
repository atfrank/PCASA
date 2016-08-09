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

#ifndef ANALYZE_H
#define ANALYZE_H

#include "Enum.hpp"
#include <vector>
#include <string>

//Forward Declaration
class Molecule;
class Coor;
class DTree;

//Abstract base class (cannot create instance of it!)
class Analyze {
  private:
    //Since this is an abstract base class
    //all members need to be accessed via a function or passed
    //directly to the analysis function or set the members as protected
    std::vector<std::string> sel;
    std::vector<Molecule*> mol;
    std::vector<double> tdata; //Time dependent data, maybe for averaging
    std::vector<std::vector<double> > fdata; //Frame data, cleared after each frame
    int ndata; //Total number of datapoints
    bool resel; //Re-do selection for each analysis, not implemented yet
    std::string ifile;
    std::string ofile;
    bool verbose;

  public:
    Analyze ();
    virtual ~Analyze ();
    void addSel(const std::string& selin);
    std::string getSel(const int& element);
    unsigned int getNSel();
    void addMol(Molecule* molin);
    void setMol(const int& element, Molecule* molin);
    void clearMol();
    void resizeNMol(const int sizein);
    Molecule* getMol(const int& element);
    unsigned int getNMol();
    void setNData(const int& ndatain);
    int& getNData();
    std::vector<double>& getTDataVec();
    std::vector<std::vector<double> >& getFDataVec();
    void setInput(const std::string& fin);
    std::string getInput();
    void setOutput (const std::string& fin);
    std::string getOutput();
    void setVerbose(bool verbosein);
    bool getVerbose();
        
    //Virtual functions
    virtual void readTopology(Molecule* molin, std::string topin="");
    virtual void setupMolSel(Molecule* molin);
    virtual void preAnalysis(Molecule* molin, std::string topin=""); 
    virtual void preAnalysis();
    virtual void runAnalysis() =0; //Pure virtual function
    virtual void postAnalysis();

    //Analysis functions
    static Coor centerOfGeometry(Molecule* mol, bool selFlag=true);
    static void averageMol(Molecule* cmpmol, Molecule* refmol, int &ndataIO);
    static double distance (const Coor& u, const Coor& v);
    static double angle (const Coor& u, const Coor& v, const Coor& w);
    static double dihedral (const Coor& t, const Coor& u, const Coor& v, const Coor& w);
    static double distance (Molecule* sel1, Molecule* sel2, bool selFlag=true);
    static double angle (Molecule* sel1, Molecule* sel2, Molecule* sel3, bool selFlag=true);
    static double dihedral (Molecule* sel1, Molecule* sel2, Molecule* sel3, Molecule* sel4,bool selFlag=true);
    static void pairwiseDistance(Molecule *mol, std::vector<std::vector<double> >& pdin);
    static void allAnglesDihedrals(Molecule *mol, std::vector<std::vector<double> >& anglesin);
    static void pcasso(Molecule* mol, std::vector<std::vector<double> > &fdataIO);
    //static std::vector<double> gyration(Molecule* mol);
};

//Derived classes

class AnalyzeCOG: public Analyze {
  public:
    void runAnalysis();
};

class AnalyzeAverage: public Analyze {
  public:
    void preAnalysis(Molecule* molin, std::string topin="");
    void runAnalysis();
    void postAnalysis();
};

class AnalyzeDistance: public Analyze {
  public:
    void setupMolSel(Molecule* molin);
    void runAnalysis();
};

class AnalyzeAngle: public Analyze {
  public:
    void setupMolSel(Molecule* molin);
    void runAnalysis();
};

class AnalyzeDihedral: public Analyze {
  private:
    DihedralEnum theta;

  public:
    AnalyzeDihedral(DihedralEnum thetain=GENERAL);
    void setDihedralType(DihedralEnum thetain=GENERAL);
    DihedralEnum getDihedralType();
    void setupMolSel(Molecule* molin);
    void runAnalysis();
};

class AnalyzePairwiseDistance: public Analyze {
  public:
    void runAnalysis();
};

class AnalyzePcasso: public Analyze {
  private:
    PcassoOutEnum pout;
    std::vector<DTree *> t;

  public:
    AnalyzePcasso(std::string delim=":");
    void setOutType(PcassoOutEnum pin);
    PcassoOutEnum getOutType();
    void preAnalysis(Molecule* molin, std::string fin="");
    void runAnalysis();
    void postAnalysis();
};

class AnalyzeCasa: public Analyze {
  private:
    std::vector<DTree *> t1;
  public:
    AnalyzeCasa(std::string delim=":");
    void preAnalysis(Molecule* molin);
    void runAnalysis();
    void runAnalysisCasa(unsigned int frame=1, std::string identification="test",bool header=false, bool predict=false, std::string fparmfile="", std::string freffile="", std::string fsasa="", const double cutoff=10.0, const int beta=-1);
};


#endif
