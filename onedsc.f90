program onedscalar
! Solve 1-D scalar advection equation
!  \partial_t \phi + U \partial_x \phi = 0


implicit none

integer             :: i,j,k
integer, parameter  :: Nx=33
real                :: x(Nx), phi(Nx), phi_n(Nx), exact(Nx)
real                :: phi_e, phi_w
real                :: dx, dt, U, Lx, Tf, t, theta, bc 
real                :: phi_nvd, beta, gamma, eps, err


Lx = 1.0
dx = Lx/(real(Nx -1))
Tf = 0.5    ! final time
U = 1.0 ! advection velocity
dt = 0.001 ! time step
bc = 1.0 ! boundary condition
beta = 0.1 ! lower limit
print*, 'CFL = ',U*dt/dx
eps = dx/100.0

!Output files
open(unit = 10, file = 'bcd.txt',status='unknown')
!open(unit = 11, file = 'exact.txt',status='unknown')
! Equally spaced mesh
do i=1,Nx
    x(i) = (i-1)*dx
enddo


!IC
t = 0.0
phi = 0.0
phi_n = 0.0
phi(1) = bc
phi_n(1) = bc


do while (t .le. Tf)
    !BC
    phi(1) = bc
    phi_n(1) = bc
    do i = 2,Nx
        ! First order upwind, U > 0
!         phi_e = phi_n(i)
!         phi_w = phi_n(i-1)

        !QUICK/2nd order upwind/Central
        ! theta = 1/8 QUICK, theta = 0 2nd order upwind, theta = 1 central (for uniform mesh only)
!         theta = 1.0/8.0
!         if (i .ne. Nx) then
!             phi_e = theta*(0.5*phi_n(i) + 0.5*phi_n(i+1)) + (1.0-theta)*(1.5*phi_n(i) - 0.5*phi_n(i-1))
!         else
!             phi_e = theta*(0.5*phi_n(i) + 0.5*phi_n(i)) + (1.0-theta)*(1.5*phi_n(i) - 0.5*phi_n(i-1))
!         endif
!         
!         if (i.ne. 2) then
!             phi_w = theta*(0.5*phi_n(i-1) + 0.5*phi_n(i)) + (1.0-theta)*(1.5*phi_n(i-1) - 0.5*phi_n(i-2))
!         else
!             phi_w = theta*(0.5*phi_n(i-1) + 0.5*phi_n(i)) + (1.0-theta)*(1.5*phi_n(i-1) - 0.5*bc)
!         endif
        
        
        !Bounded central difference (NVD)
        !face e
        phi_nvd = (phi_n(i) - phi_n(i-1))/(phi_n(i+1) - phi_n(i-1) + eps)
        ! find blending function
        gamma = phi_nvd/beta
        
        if ((phi_nvd .gt. 1.0) .or. (phi_nvd .lt. 0.0)) then !not bounded, use first order upwind
            phi_e = phi_n(i)
        else if ((phi_nvd .ge. beta) .and. (phi_nvd .lt. 1.0)) then ! use central scheme
            if (i .eq. Nx) then
                phi_e = phi_n(i)
            else
                phi_e = 0.5*(phi_n(i) + phi_n(i+1))
            endif
        else ! use blended CD/UD
            if (i .eq. Nx) then
                phi_e = phi_n(i)
            else
                phi_e = (1.0 - gamma*0.5)*phi_n(i) + gamma*0.5*phi_n(i+1)
            endif
        endif
        
        !face w
        phi_nvd = (phi_n(i-1) - phi_n(i-2))/(phi_n(i) - phi_n(i-2) + eps)
        ! find blending function
        gamma = phi_nvd/beta
        
        if ((phi_nvd .gt. 1.0) .or. (phi_nvd .lt. 0.0)) then !not bounded, use first order upwind
            phi_w = phi_n(i-1)
        else if ((phi_nvd .ge. beta) .and. (phi_nvd .lt. 1.0)) then ! use central scheme
            phi_w = 0.5*(phi_n(i-1) + phi_n(i))
        else ! use blended CD/UD
            phi_w = (1.0 - gamma*0.5)*phi_n(i-1) + gamma*0.5*phi_n(i)
        endif
        
        
        !Explicit update
        phi(i) = phi_n(i) - U*dt/dx * (phi_e - phi_w)
        phi_n(i) = phi(i)
    enddo
    t = t + dt
enddo


err = 0.0
do i=1,Nx
    if (x(i) .le. U*t) then
        exact(i) = U
    else
        exact(i) = 0.0
    endif
    write(10,*) x(i), phi(i)
    !write(11,*) x(i), exact(i)
    
    err = err + (phi(i) - exact(i))**2.0
enddo
err = sqrt(err/Nx)

print*, err
close(10)
!close(11)
end program onedscalar
