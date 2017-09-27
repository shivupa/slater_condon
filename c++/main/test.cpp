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

#include <cmath>
#include <chrono>
#include <vector>
#include <iostream>
#include <bitset>


#include "../read_input.h"
#include "../rdtsc.c"
#include "../nsubst1.x.h"
#include "../exc.h"

int main(int argc, char **argv){
    const int Nint = 2;
    const int ndet = 10000;
    const int mo_num = 105;

    std::vector<double> density_matrix; //= new double[mo_num*mo_num];
    std::vector<double> coef; //= new double[ndet];
    std::vector< std::vector< std::vector< unsigned int > > > det; // = new std::vector<unsigned long long int>[ndet*2*Nint];
    int n_excit = 0;

    read_input(det, coef, Nint, ndet, mo_num,  "h2o_determinants.dat", "cu.coef");
    std::cout << std::endl;
    std::cout << "ndet = " << ndet << std::endl;
    std::cout << std::endl;

    //-------------------------------------------------
    double t = 0.0;
    double t0 = 0.0;
    double t1 = 0.0;
    double t2 = 0.0;
    double error = 0.0;

    std::chrono::high_resolution_clock::time_point cpu0=std::chrono::high_resolution_clock::now();
    double events = 0.0;
    int res;
    for(int l = 0; l < ndet; l++){
        t0 = irp_rdtsc_();
        for(int k = 0; k < ndet; k++){
            res = n_excitations(det, l, k, Nint);
        }
        t1 = irp_rdtsc_();
        t += (t1-t0);
        t2 *= (t1-t0)*(t1-t0);
    }
    std::cout << "test" << std::endl;

    events = ndet*ndet;
    std::chrono::high_resolution_clock::time_point cpu1=std::chrono::high_resolution_clock::now();
    error = std::sqrt(std::abs(std::pow((t/events),2)-t2/events)/events );

    std::cout << "Cycles n_excitations : " << t/events << " +/- " << error/std::sqrt(events) << std::endl;
    std::cout << "CPU    n_excitations : " << std::chrono::duration_cast<std::chrono::microseconds>(cpu1-cpu0).count()/1e6 << " s"<< std::endl;
    std::cout << res << std::endl;
    //-------------------------------------------------
    //-------------------------------------------------
    t = 0.0;
    t2 = 0.0;
    cpu0=std::chrono::high_resolution_clock::now();
    std::vector< std::vector< std::vector<unsigned int > > > exc (3,std::vector< std::vector< unsigned int > >(2,std::vector <unsigned int > (0,0)));
    int degree;
    double phase = 0.0;
    for(int l = 0; l < ndet; l++){
        t0 = irp_rdtsc_();
        for(int k = 0; k < ndet; k++){
            //get excitation
            get_excitation(det, phase, exc, degree, l, k, Nint);
        }
    }
}
/*
    integer    :: i,k,l
    double precision :: phase
    integer    :: exc(0:2,2,2)
    double precision :: t0, t1, irp_rdtsc
    double precision :: t, t2, t3(-1:2), t4(-1:2)
    double precision :: nb(-1:2)
    double precision :: cpu0, cpu1
    double precision :: events
    double precision :: error
    integer,parameter  :: lmax = 10000
    integer:: res

    //-------------------------------------------------

    t=0.d0
    t3=0.d0
    t4=0.d0
    do l=1,ndet-1

     t0 = irp_rdtsc()
     do k=1,lmax
      //DIR$ FORCEINLINE
      call get_excitation(det(1,1,l),det(1,1,1),exc,i,phase,N_int)
     enddo
     if (i==1) then
     do k=lmax+1,lmax*100
      //DIR$ FORCEINLINE
      call get_excitation(det(1,1,l),det(1,1,1),exc,i,phase,N_int)
     enddo
     endif
     t1 = irp_rdtsc()
     t3(i) = t3(i)+(t1-t0)
     t4(i) = t4(i)+(t1-t0)*(t1-t0)
     nb(i) = nb(i) + (k-1.d0)

     t0 = irp_rdtsc()
     do k=1,lmax
      //DIR$ FORCEINLINE
      call get_excitation(det(1,1,l),det(1,1,l+1),exc,i,phase,N_int)
     enddo
     if (i==1) then
     do k=lmax+1,lmax*100
      //DIR$ FORCEINLINE
      call get_excitation(det(1,1,l),det(1,1,l+1),exc,i,phase,N_int)
     enddo
     endif
     t1 = irp_rdtsc()
     t3(i) = t3(i)+(t1-t0)
     t4(i) = t4(i)+(t1-t0)*(t1-t0)
     nb(i) = nb(i) + (k-1.d0)


     t0 = irp_rdtsc()
     do k=1,lmax
      //DIR$ FORCEINLINE
      call get_excitation(det(1,1,l),det(1,1,ndet-l+1),exc,i,phase,N_int)
     enddo
     if (i==1) then
     do k=lmax+1,lmax*100
      //DIR$ FORCEINLINE
      call get_excitation(det(1,1,l),det(1,1,ndet-l+1),exc,i,phase,N_int)
     enddo
     endif
     t1 = irp_rdtsc()
     t3(i) = t3(i)+(t1-t0)
     t4(i) = t4(i)+(t1-t0)*(t1-t0)
     nb(i) = nb(i) + (k-1.d0)


     t0 = irp_rdtsc()
     do k=1,lmax/10
      //DIR$ FORCEINLINE
      call get_excitation(det(1,1,l),det(1,1,l),exc,i,phase,N_int)
     enddo
     t1 = irp_rdtsc()
     t3(0) = t3(0)+(t1-t0)
     t4(0) = t4(0)+(t1-t0)*(t1-t0)
     nb(0) = nb(0) + (k-1.d0)
    enddo
    error = sqrt( abs((t3(0)/nb(0))**2-t4(0)/nb(0))/nb(0) )
    print *,  'Cycles zero   :',  t3(0)/nb(0), ' +/- ', error/sqrt(nb(0))//, nb(0)
    error = sqrt( abs((t3(1)/nb(1))**2-t4(1)/nb(1))/nb(1) )
    print *,  'Cycles single :',  t3(1)/nb(1), ' +/- ', error/sqrt(nb(1))//, nb(1)
    error = sqrt( abs((t3(2)/nb(2))**2-t4(2)/nb(2))/nb(2) )
    print *,  'Cycles double :',  t3(2)/nb(2), ' +/- ', error/sqrt(nb(2))//, nb(2)
    error = sqrt( abs((t3(-1)/nb(-1))**2-t4(-1)/nb(-1))/nb(-1) )
    print *,  'Cycles other  :',  t3(-1)/nb(-1), ' +/- ', error/sqrt(nb(-1))//, nb(-1)
    //-------------------------------------------------


    call cpu_time(cpu0)
    t0 = irp_rdtsc()
    call compute_density_matrix(det,ndet,coef,mo_num, &
                   N_int,density_matrix)
    t1 = irp_rdtsc()
    call cpu_time(cpu1)
    print *,  'Cycles density matrix : ', (t1-t0)/(ndet*(ndet-1)/2)
    print *,  'CPU    density matrix : ', cpu1-cpu0
   //print *, (density_matrix(k,k), k=1,mo_num)

  end //----------

subroutine print_key( key, N_int  )
//////// print the physical meaning of "key(N_int,2)" in an explicit way
 implicit none
 integer :: N_int
 integer*8, intent(in)    :: key(N_int,2)
 integer, parameter       :: alpha = 1
 integer, parameter       :: beta  = 2

 integer :: i, j, k, ibuf
 integer*8 :: itemp
 character*(1) :: buffer(10000)

 ibuf = 1
 do k=alpha,beta
   do i=1,N_int
    itemp = 1
    do j=1,64
      if (iand(itemp,key(i,k)) == itemp) then
        buffer(ibuf) = '1'
      else
        buffer(ibuf) = '.'
      endif
      ibuf = ibuf+1
      itemp = ishft(itemp,1)
    enddo
   enddo
 enddo
 print *, key(:,1), key(:,2)
 print *, buffer(1:ibuf)
 print *,  ''

end
*/
