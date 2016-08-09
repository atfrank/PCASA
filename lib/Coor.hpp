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

#ifndef VECTOR_H
#define VECTOR_H

//#include <cmath>
//#include <iostream>

class Coor {
  private:
    double xcoor;
    double ycoor;
    double zcoor;

  public:
    Coor();
    Coor(double xcoorin, double ycoorin, double zcoorin); //Constructor
    Coor(const Coor& vec); //Overload Constructor 

    Coor& operator= (const Coor& vec);
    Coor& operator= (const double val);
    //Addition
    Coor operator+ (const Coor& vec) const;
    Coor& operator+= (const Coor& vec);
    Coor operator+ (const double val) const;
    Coor& operator+= (const double val);
    //Subtraction
    Coor operator- (const Coor& vec) const;
    Coor& operator-= (const Coor& vec);
    Coor operator- (const double val) const;
    Coor& operator-= (const double val);
    //Multiplication
    Coor operator* (const Coor& vec) const;
    Coor& operator*= (const Coor& vec);
    Coor operator* (const double val) const;
    Coor& operator*= (const double val);
    //Division
    Coor operator/ (const Coor& vec) const;
    Coor& operator/= (const Coor& vec);
    Coor operator/ (const double val) const;
    Coor& operator/= (const double val);

    double& x(){return xcoor;};
    double& y(){return ycoor;};
    double& z(){return zcoor;};

    Coor operator- () const;
    double dot (const Coor& vec) const; //Dot Product
    Coor cross (const Coor& vec) const; //Cross Product
    double norm () const;
};

#endif
