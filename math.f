c***********************************************************************
c***********************************************************************
      subroutine evolvepatch(
     &phi_p , pi_p , phi_c, pi_c,
     &phigi, phigj, phigk,
     &pigi, pigj, pigk ,
     &ifirst, ilast, jfirst, jlast , kfirst, klast,
     &dx, dt )

      implicit none
      integer ifirst, ilast, jfirst, jlast , kfirst, klast
      integer pigi, pigj, pigk ,
     &        phigi, phigj, phigk
      double precision phi_p(ifirst-phigi:ilast+phigi,
     &                     jfirst-phigj:jlast+phigj,
     &                     kfirst-phigk:jlast+phigk)
      double precision pi_p(ifirst-pigi:ilast+pigi,
     &                     jfirst-pigj:jlast+pigj,
     &                     kfirst-pigk:jlast+pigk)
      double precision phi_c(ifirst-phigi:ilast+phigi,
     &                     jfirst-phigj:jlast+phigj,
     &                     kfirst-phigk:jlast+phigk)
      double precision pi_c(ifirst-pigi:ilast+pigi,
     &                     jfirst-pigj:jlast+pigj,
     &                     kfirst-pigk:jlast+pigk)
      double precision dx(0:2)
      double precision dxi, dyi, dzi, dxi2, dyi2, dzi2
      double precision dt, lap
      integer i, j, k
      

      
      dxi = 1./dx(0)
      dyi = 1./dx(1)
      dzi = 1./dx(2)
      dxi2 = dxi*dxi
      dyi2 = dyi*dyi
      dzi2 = dzi*dzi

      do k=kfirst,klast
         do j=jfirst,jlast
            do i=ifirst,ilast
               lap =
     &         dxi2*(phi_p(i+1,j,k)+phi_p(i-1,j,k) - 2.0 * phi_p(i,j,k))
     &       + dyi2*(phi_p(i,j+1,k)+phi_p(i,j-1,k) - 2.0 * phi_p(i,j,k))
     &       + dzi2*(phi_p(i,j,k+1)+phi_p(i,j,k-1) - 2.0 * phi_p(i,j,k))
               phi_c(i,j,k) = phi_p(i,j,k) + pi_p(i,j,k)*dt
               pi_c(i,j,k) = pi_p(i,j,k) + lap*dt
            enddo
         enddo
      enddo
      
      return
      end
c***********************************************************************
