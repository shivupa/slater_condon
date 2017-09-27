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
#include <iostream>

int n_excitations(std::vector< std::vector< std::vector< unsigned int > > >& det, int det1, int det2, int Nint){
    int excitations = 0;

    for(int l = 0; l < Nint; l++){
        // excitations += (int)(det[det1][l][0] ^ det[det2][l][0]).count(); if using bitsets
        // excitations += (int)(det[det1][l][1] ^ det[det2][l][1]).count(); if using bitsets
        excitations += __builtin_popcount((det[det1][l][0] ^ det[det2][l][0]));
        excitations += __builtin_popcount((det[det1][l][1] ^ det[det2][l][1]));
        // excitations += __builtin_popcountl((det[det1][l][0] ^ det[det2][l][0]).to_ulong());
        // excitations += __builtin_popcountl((det[det1][l][1] ^ det[det2][l][1]).to_ulong());
        // I'm not sure if this will be faster. I think they should both end up calling popcnt
        // intrinsics but i havent checked.

    }
    return excitations<<1;
}
