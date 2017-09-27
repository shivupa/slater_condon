//   Slater-Condon Rules
//   Copyright (C) 2013 Anthony Scemama <scemama@irsamc.ups-tlse.fr>
//                      Emmanuel Giner <emmanuel_giner_jr@hotmail.fr>
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License along
//   with this program; if not, write to the Free Software Foundation, Inc.,
//   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

void read_input(std::vector< std::vector< std::vector< unsigned int> > >& det,  std::vector<double>& coef, const int Nint, const int ndet, const int mo_num, std::string filename, std::string filename2){
    std::ifstream infile(filename);
    int Nint_from_file;
    infile >> Nint_from_file;
    assert(Nint_from_file == Nint);
    unsigned int val;
    std::vector< unsigned int  > spin;
    std::vector< std::vector< unsigned int > > one_det;

    for( int i = 0 ; i < ndet; i++){
        one_det.clear();
        for( int j = 0; j < Nint; j++){
            spin.clear();
            for( int k = 0; k < 2; k++){
                infile >> val;
                spin.push_back(val);
            }
            one_det.push_back(spin);
        }
        det.push_back(one_det);
    }
    std::cout << "Shape dets : (" << det.size() << "," << det[0].size() << "," << det[0][0].size() << ")" << std::endl;
    // while (infile >> val){
    //     det.push_back(val);
    // }
    // infile.close();

    std::ifstream infile2(filename);
    double val2;
    while (infile2 >> val2){
        coef.push_back(val2);
    }
    infile2.close();
}
