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
#include <algorithm>
#include <vector>

void get_excitation(std::vector< std::vector< std::vector< unsigned int > > >& det, double& phase, std::vector< std::vector< std::vector< unsigned int > > >& exc, int& degree, int det1, int det2, int Nint);

void get_double_excitation(std::vector< std::vector< std::vector< unsigned int > > >& det, double& phase, std::vector< std::vector< std::vector< unsigned int > > >& exc, int det1, int det2, int Nint);

void get_single_excitation(std::vector< std::vector< std::vector<  unsigned int > > >& det, double& phase, std::vector< std::vector< std::vector< unsigned int > > >& exc, int det1, int det2, int Nint);
