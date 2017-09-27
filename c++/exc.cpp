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
#include "nsubst1.x.h"
#include "exc.h"

void get_excitation(std::vector< std::vector< std::vector< unsigned int > > >& det, double& phase, std::vector< std::vector< std::vector< unsigned int > > >& exc, int& degree, int det1, int det2, int Nint){

    int n_excit = n_excitations(det, det1, det2, Nint);

    switch (n_excit){
        case 3 :
            degree = -1;
            return;
        case 2 :
            get_double_excitation(det, phase, exc, det1, det2, Nint);
            return;
        case 1 :
            get_single_excitation(det, phase, exc, det1, det2, Nint);
            return;
        default :
            return;
    }
}
/*
subroutine get_excitation(det1,det2,exc,degree,phase,Nint)
 implicit none
 integer, intent(in)  :: Nint
 integer*8, intent(in)  :: det1(Nint,2), det2(Nint,2)
 integer, intent(out) :: exc(0:2,2,2)
 integer, intent(out) :: degree
 double precision, intent(out) :: phase

 integer :: n_excitations

 degree = n_excitations(det1,det2,Nint)

 select case (degree)

   case (3:)
     degree = -1
     return

   case (2)
     call get_double_excitation(det1,det2,exc,phase,Nint)
     return

   case (1)
     call get_single_excitation(det1,det2,exc,phase,Nint)
     return

   case(0)
     return

   end select
end
*/
void get_double_excitation(std::vector< std::vector< std::vector< unsigned int > > >& det, double& phase, std::vector< std::vector< std::vector< unsigned int > > >& exc, int det1, int det2, int Nint){
    double* phase_dble[2];
    *phase_dble[0] = -1;
    *phase_dble[1] = 1;

    exc[0][0][0] = 0;
    exc[0][1][0] = 0;
    exc[0][0][1] = 0;
    exc[0][1][1] = 0;

    unsigned int particle;
    unsigned int hole;
    int ishift;
    unsigned int tz;
    unsigned int tmp;

    int idx_particle;
    int idx_hole;
    int nexc = 0;
    int nperm = 0;
    for(int ispin = 0 ; ispin < 2; ispin++){
        idx_particle = 0;
        idx_hole = 0;
        ishift = -63;
        for(int l = 0; l < Nint; l++){
            ishift += 64;
            if(det[det1][l][ispin] != det[det2][l][ispin]){
                tmp = det[det1][l][ispin] ^ det[det2][l][ispin];
                particle = tmp & det[det2][l][ispin];
                hole = tmp & det[det1][l][ispin];
                while (particle != 0){
                    tz = __builtin_ctz(particle);
                    nexc += 1;
                    exc[0][1][ispin] += 1;
                    exc[idx_particle][1][ispin] = tz+ishift ;
                    idx_particle += 1;
                    particle = particle & (particle-1);
                }
                while (hole != 0){
                    tz = __builtin_ctz(hole);
                    nexc += 1;
                    exc[0][0][ispin] += 1;
                    exc[idx_hole][0][ispin] = tz+ishift ;
                    idx_hole +=1;
                    hole = hole & (hole-1);
                }
            }
            if(nexc == 4){
                return;
            }
        }
        for (int i = 0; i < exc[0][1][ispin]; i++){
            unsigned int low  = std::min(exc[i][0][ispin], exc[i][1][ispin]);
            unsigned int high = std::max(exc[i][0][ispin], exc[i][1][ispin]);
            unsigned int j = ((low-1) >> 6) + 1;
            unsigned int n = (low-1) & 63;
            unsigned int k = ((high-1) >> 6) +1;
            unsigned int m = (high-1) & 63;
            if (j == k){
                nperm += __builtin_popcount((det[det1][j][ispin])  & ((~(1 << (n+1)))+1 & (1 << m) - 1)) ;
            } else {
                nperm = __builtin_popcount((det[det1][k][ispin]) & ((1 << m) - 1)) + __builtin_popcount((det[det1][j][ispin]) & ((~(1 << (n+1)))+1));
                for(int i = (j+1); i<(k-1); i++){
                    nperm += __builtin_popcount(det[det1][i][ispin]);
                }
            }
        }
        if (exc[0][1][ispin] == 2){
            unsigned int a = std::min(exc[1][0][ispin], exc[1][1][ispin]);
            unsigned int b = std::max(exc[1][0][ispin], exc[1][1][ispin]);
            unsigned int c = std::min(exc[2][0][ispin], exc[2][1][ispin]);
            unsigned int d = std::max(exc[2][0][ispin], exc[2][1][ispin]);
            if ((c>a) && (c<b) && (d>b)){
                nperm+=1;
                //exit()
            }
        }
        phase = *phase_dble[(int)(nperm & 1)];
    }
}
/*

subroutine get_double_excitation(det1,det2,exc,phase,Nint)
 implicit none
 integer, intent(in)  :: Nint
 integer*8, intent(in)  :: det1(Nint,2), det2(Nint,2)
 integer, intent(out) :: exc(0:2,2,2)
 double precision, intent(out) :: phase
 integer :: l, ispin, idx_hole, idx_particle, ishift
 integer :: i,j,k,m,n,high, low,a,b,c,d,nperm,tz,nexc
 integer*8 :: hole, particle, tmp
 double precision, parameter :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)
 exc(0,1,1) = 0
 exc(0,2,1) = 0
 exc(0,1,2) = 0
 exc(0,2,2) = 0
 nexc=0
 nperm=0
 do ispin = 1,2
  idx_particle = 0
  idx_hole = 0
  ishift = -63
  do l=1,Nint
   ishift = ishift + 64
   if (det1(l,ispin) == det2(l,ispin))  then
     cycle
   end if
   tmp = xor( det1(l,ispin), det2(l,ispin) )
   particle = iand(tmp, det2(l,ispin))
   hole     = iand(tmp, det1(l,ispin))
   do while (particle /= 0_8)
     tz = trailz(particle)
     nexc = nexc+1
     idx_particle = idx_particle + 1
     exc(0,2,ispin) = exc(0,2,ispin) + 1
     exc(idx_particle,2,ispin) = tz+ishift
     particle = iand(particle,particle-1_8)
   end do
   do while (hole /= 0_8)
     tz = trailz(hole)
     nexc = nexc+1
     idx_hole = idx_hole + 1
     exc(0,1,ispin) = exc(0,1,ispin) + 1
     exc(idx_hole,1,ispin) = tz+ishift
     hole = iand(hole,hole-1_8)
   end do
   if (nexc == 4) exit
  end do

  do i=1,exc(0,1,ispin)
    low  = min(exc(i,1,ispin),exc(i,2,ispin))
    high = max(exc(i,1,ispin),exc(i,2,ispin))
    j = ishft(low-1,-6)+1
    n = iand(low-1,63)
    k = ishft(high-1,-6)+1
    m = iand(high-1,63)
    if (j==k) then
      nperm = nperm + popcnt(iand(det1(j,ispin),  &
         iand( not(ishft(1_8,n+1))+1 ,ishft(1_8,m)-1)))
    else
      nperm = nperm + popcnt(iand(det1(k,ispin),  &
                             ishft(1_8,m)-1)) &
                    + popcnt(iand(det1(j,ispin),  &
                             not(ishft(1_8,n+1))+1))
      do l=j+1,k-1
        nperm = nperm + popcnt(det1(l,ispin))
      end do
    end if
  end do
  if (exc(0,1,ispin) == 2) then
    a = min(exc(1,1,ispin), exc(1,2,ispin))
    b = max(exc(1,1,ispin), exc(1,2,ispin))
    c = min(exc(2,1,ispin), exc(2,2,ispin))
    d = max(exc(2,1,ispin), exc(2,2,ispin))
    if (c>a .and. c<b .and. d>b) nperm = nperm + 1
    exit
   end if
 end do
 phase = phase_dble(iand(nperm,1))

end

*/
void get_single_excitation(std::vector< std::vector< std::vector<  unsigned int > > >& det, double& phase, std::vector< std::vector< std::vector< unsigned int > > >& exc, int det1, int det2, int Nint){
    double* phase_dble[2];
    *phase_dble[0] = -1;
    *phase_dble[1] = 1;

    exc[0][0][0] = 0;
    exc[0][1][0] = 0;
    exc[0][0][1] = 0;
    exc[0][1][1] = 0;

    int particle;
    int hole;
    int ishift;
    int tz;
    int tmp;
    int nperm;

    for(int ispin = 0 ; ispin < 2; ispin++){
        ishift = -63;
        for(int l = 0; l < Nint; l++){
            ishift += 64;
            if(det[det1][l][ispin] != det[det2][l][ispin]){
                tmp = det[det1][l][ispin] ^ det[det2][l][ispin];
                particle = tmp & det[det2][l][ispin];
                hole = tmp & det[det1][l][ispin];
                if (particle != 0){
                    tz = __builtin_ctz(particle);
                    exc[0][1][ispin] = 1;
                    exc[1][1][ispin] = tz+ishift;
                }
                if (hole != 0){
                    tz = __builtin_ctz(hole);
                    exc[0][0][ispin] = 1;
                    exc[1][0][ispin] = tz+ishift ;
                }
                if((exc[0][0][ispin] & exc[0][1][ispin]) == 1){
                    unsigned int low  = std::min(exc[0][0][ispin], exc[0][1][ispin]);
                    unsigned int high = std::max(exc[0][0][ispin], exc[0][1][ispin]);
                    unsigned int j = ((low-1) >> 6) + 1;
                    unsigned int n = (low-1) & 63;
                    unsigned int k = ((high-1) >> 6) +1;
                    unsigned int m = (high-1) & 63;
                    if (j == k){
                        nperm = __builtin_popcount((det[det1][j][ispin])  & ((~(1 << (n+1)))+1 & (1 << m) - 1)) ;
                    } else {
                        nperm = __builtin_popcount((det[det1][k][ispin]) & ((1 << m) - 1)) + __builtin_popcount((det[det1][j][ispin]) & ((~(1 << (n+1)))+1));
                        for(int i = (j+1); i<(k-1); i++){
                            nperm += __builtin_popcount(det[det1][i][ispin]);
                        }
                    }
                    phase = *phase_dble[nperm & 1];
                }
            }
        }
    }
}
/*
subroutine get_single_excitation(det1,det2,exc,phase,Nint)
 implicit none
 integer, intent(in)  :: Nint
 integer*8, intent(in)  :: det1(Nint,2)
 integer*8, intent(in)  :: det2(Nint,2)
 integer, intent(out) :: exc(0:2,2,2)
 double precision, intent(out) :: phase
 integer :: tz, l, ispin, ishift, nperm, i, j, k, m, n, high, low
 integer*8 :: hole, particle, tmp
 double precision, parameter :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)

 exc(0,1,1) = 0
 exc(0,2,1) = 0
 exc(0,1,2) = 0
 exc(0,2,2) = 0
 do ispin = 1,2
   ishift = -63
   do l=1,Nint
     ishift = ishift + 64
     if (det1(l,ispin) == det2(l,ispin)) cycle
     tmp = xor( det1(l,ispin), det2(l,ispin) )
     particle = iand(tmp, det2(l,ispin))
     hole     = iand(tmp, det1(l,ispin))
     if (particle /= 0_8) then
       tz = trailz(particle)
       exc(0,2,ispin) = 1
       exc(1,2,ispin) = tz+ishift
     end if
     if (hole /= 0_8) then
       tz = trailz(hole)
       exc(0,1,ispin) = 1
       exc(1,1,ispin) = tz+ishift
     end if

     if ( iand(exc(0,1,ispin),exc(0,2,ispin)) == 1 ) then
       low  = min(exc(1,1,ispin),exc(1,2,ispin))
       high = max(exc(1,1,ispin),exc(1,2,ispin))
       j = ishft(low-1,-6)+1
       n = iand(low-1,63)
       k = ishft(high-1,-6)+1
       m = iand(high-1,63)
       if (j==k) then
         nperm = popcnt(iand(det1(j,ispin), &
            iand( not(ishft(1_8,n+1))+1_8 ,ishft(1_8,m)-1_8)))
       else
         nperm = popcnt(iand(det1(k,ispin), ishft(1_8,m)-1_8)) + &
                 popcnt(iand(det1(j,ispin), not(ishft(1_8,n+1))+1_8))
         do i=j+1,k-1
           nperm = nperm + popcnt(det1(i,ispin))
         end do
       end if
       phase = phase_dble(iand(nperm,1))
       return
     end if
   end do
 end do
end
*/
