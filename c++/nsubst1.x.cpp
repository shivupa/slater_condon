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
#include <vector>
#include <bitset>
#include <iostream>

int n_excitations(std::vector< std::vector< std::vector< std::bitset<64> > > >& det, int det1, int det2, int Nint){
    int excitations = 0;

    for(int l = 0; l < Nint; l++){
        excitations += (int)(det[det1][0][l] ^ det[det2][0][l]).to_ulong();
        excitations += (int)(det[det1][1][l] ^ det[det2][1][l]).to_ulong();
        // excitations += __builtin_popcountl((det[det1][0][l] ^ det[det2][0][l]).to_ulong());
        // excitations += __builtin_popcountl((det[det1][1][l] ^ det[det2][1][l]).to_ulong());

    }
    return excitations<<=1;
}
