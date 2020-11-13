module constant
    implicit none
    double precision,parameter :: A=0.10
    double precision,parameter :: B=0.28
    double precision,parameter :: C=0.015
    double precision,parameter :: D=0.06
    double precision,parameter :: E0=0.05 
    double precision,parameter :: dt=1.0D-1
    double precision,parameter :: m=2000
    double precision,parameter :: pi=3.1415926
end module constant

!********************************************************
!              Preparation for dynamics                 *
!********************************************************
module Hamiltonian
    use constant
    implicit none
contains
!********************************************************
!                  Diabatic Hamiltonian                 *
!********************************************************
    double precision function V11(x)
        implicit none
        double precision :: x
        v11=0
        return
    endfunction

    double precision function V22(x)
        implicit none
        double precision :: x
        V22=-A*exp(-B*x**2)+E0
        return
    endfunction

    double precision function V12(x)
        implicit none
        double precision :: x
        V12=C*exp(-D*x**2)
        return 
    endfunction


!********************************************************
!          derivative of Diabatic Hamiltonian           *
!********************************************************
    double precision function D_V11(x)
        implicit none
        double precision :: x
        D_V11=0
        return
    endfunction

    double precision function D_V22(x)
        implicit none
        double precision :: x
        D_V22=2*A*B*x*exp(-B*x**2)
        return
    endfunction

    double precision function D_V12(x)
        implicit none
        double precision :: x
        D_V12=-2*C*D*x*Exp(-D*x**2)
        return
    endfunction

!**************************************************
!            Adiabatic Hamiltonian                !
!**************************************************
    double precision function E1(x)
        implicit none
        double precision :: x
        E1=0.5*(V11(x)+V22(x)-sqrt((V11(x)-V22(x))**2+4*V12(x)**2))
        return
    endfunction

    double precision function E2(x)
        implicit none
        double precision :: x
        E2=0.5*(V11(x)+V22(x)+sqrt((V11(x)-V22(x))**2+4*V12(x)**2))
        return
    endfunction

!**************************************************
!      derivative of Adiabatic Hamiltonian        !
!**************************************************
double precision function D_E1(x)
    implicit none
    double precision :: x
    D_E1=0.5*(D_V11(x)+D_V22(x)-((V11(x)-V22(x))*(D_V11(x)-D_V22(x))+4*V12(x)*D_V12(x))/sqrt((V11(x)-V22(x))**2+4*V12(x)**2))
    return
endfunction

double precision function D_E2(x)
    implicit none
    double precision :: x
    D_E2=0.5*(D_V11(x)+D_V22(x)+((V11(x)-V22(x))*(D_V11(x)-D_V22(x))+4*V12(x)*D_V12(x))/sqrt((V11(x)-V22(x))**2+4*V12(x)**2))
    return
endfunction

!**************************************************
!            Non-adiabatic coupling               !
!**************************************************
    double precision function d12(x)
        implicit none
        double precision :: x
        double precision :: c1(2),c2(2)
        double precision :: norm1,norm2

        c1(1)=-(-V11(x)+V22(x)+sqrt(V11(x)**2+4*V12(x)**2-2*V11(x)*V22(x)+V22(x)**2))/(2*V12(x))
        c1(2)=1

        c2(1)=-(-V11(x)+V22(x)-sqrt(V11(x)**2+4*V12(x)**2-2*V11(x)*V22(x)+V22(x)**2))/(2*V12(x))
        c2(2)=1

        norm1=sqrt(c1(1)**2+c1(2)**2)
        norm2=sqrt(c2(1)**2+c2(2)**2)

        c1(:)=c1(:)/norm1
        c2(:)=c2(:)/norm2

        !calculate d12
        d12=-(c1(1)*(D_V11(x)*c2(1)+D_V12(x)*c2(2))+c1(2)*(D_V12(x)*c2(1)+D_V22(x)*c2(2)))/(E2(x)-E1(x))

        return
    endfunction

    double precision function d21(x)
        implicit none
        double precision :: x
        d21=-d12(x)
        return
    endfunction

end module Hamiltonian

!********************************************************
!         Dynamics of electron and nuclear              *
!********************************************************
module dynamics
    use constant
    use Hamiltonian
    implicit none
contains

!evolution function for c1
complex(kind=8) function f(t,c1,c2,x,p,activeindex)                
    implicit none
    complex(kind=8) :: c1,c2
    double precision :: t,x,p
    integer :: activeindex
    double precision :: pna

    if (activeindex==1) then
        if ((p**2/(2*m)+E1(x)-E2(x))>0) then
            pna=sqrt(1+(2*m*(E1(x)-E2(x)))/(p**2))*p
        elseif ((p**2/(2*m)+E1(x)-E2(x))<0) then
            pna=0
        endif
        f=(0,1.)*(p**2/m)*c1-(p/m)*d12(x)*c2
    elseif (activeindex==2) then
        pna=sqrt(1+(2*m*(E2(x)-E1(x)))/(p**2))*p
        f=(0,1.)*(p*pna/m)*c1-(p/m)*d12(x)*c2
    endif

    return
endfunction

!evolution function for c2
complex(kind=8) function g(t,c1,c2,x,p,activeindex)                
    implicit none
    complex(kind=8) :: c1,c2
    double precision :: t,x,p
    integer :: activeindex
    double precision :: pna

    if (activeindex==1) then
        if ((p**2/(2*m)+E1(x)-E2(x))>0) then
            pna=sqrt(1+(2*m*(E1(x)-E2(x)))/(p**2))*p
        elseif ((p**2/(2*m)+E1(x)-E2(x))<0) then
            pna=0
        endif
        g=(0,1.)*(p*pna/m)*c2-(p/m)*d21(x)*c1
    elseif (activeindex==2) then
        pna=sqrt(1+(2*m*(E2(x)-E1(x)))/(p**2))*p
        g=(0,1.)*(p**2/m)*c2-(p/m)*d21(x)*c1
    endif

    return
endfunction

!evolution function for x
double precision function q(t,c1,c2,x,p)
    implicit none
    complex(kind=8) :: c1,c2
    double precision :: t,x,p
    q=p/m
    return
endfunction

!evolution function for p
double precision function r(t,c1,c2,x,p,activeindex)
    implicit none
    complex(kind=8) :: c1,c2
    double precision :: t,x,p
    integer :: activeindex
    
    if (activeindex==1) then
        r=-D_E1(x)
    elseif (activeindex==2) then
        r=-D_E2(x)
    endif
    return
endfunction

!********************************************************
!                       RK4                             *
!********************************************************
    subroutine RK4(t,c1,c2,x,p,activeindex)
        implicit none
        complex(kind=8) :: c1,c2
        double precision :: t,x,p
        complex(kind=8) :: K1,K2,K3,K4
        complex(kind=8) :: L1,L2,L3,L4
        double precision :: M1,M2,M3,M4
        double precision :: N1,N2,N3,N4
        integer :: activeindex

        K1=f(t,c1,c2,x,p,activeindex)
        L1=g(t,c1,c2,x,p,activeindex)
        M1=q(t,c1,c2,x,p)
        N1=r(t,c1,c2,x,p,activeindex)

        K2=f(t+(dt/2),c1+(dt/2)*K1,c2+(dt/2)*L1,x+(dt/2)*M1,p+(dt/2)*N1,activeindex)
        L2=g(t+(dt/2),c1+(dt/2)*K1,c2+(dt/2)*L1,x+(dt/2)*M1,p+(dt/2)*N1,activeindex)
        M2=q(t+(dt/2),c1+(dt/2)*K1,c2+(dt/2)*L1,x+(dt/2)*M1,p+(dt/2)*N1)
        N2=r(t+(dt/2),c1+(dt/2)*K1,c2+(dt/2)*L1,x+(dt/2)*M1,p+(dt/2)*N1,activeindex)

        K3=f(t+(dt/2),c1+(dt/2)*K2,c2+(dt/2)*L2,x+(dt/2)*M2,p+(dt/2)*N2,activeindex)
        L3=g(t+(dt/2),c1+(dt/2)*K2,c2+(dt/2)*L2,x+(dt/2)*M2,p+(dt/2)*N2,activeindex)
        M3=q(t+(dt/2),c1+(dt/2)*K2,c2+(dt/2)*L2,x+(dt/2)*M2,p+(dt/2)*N2)
        N3=r(t+(dt/2),c1+(dt/2)*K2,c2+(dt/2)*L2,x+(dt/2)*M2,p+(dt/2)*N2,activeindex)

        K4=f(t+dt,c1+dt*K3,c2+dt*L3,x+dt*M3,p+dt*N3,activeindex)
        L4=g(t+dt,c1+dt*K3,c2+dt*L3,x+dt*M3,p+dt*N3,activeindex)
        M4=q(t+dt,c1+dt*K3,c2+dt*L3,x+dt*M3,p+dt*N3)
        N4=r(t+dt,c1+dt*K3,c2+dt*L3,x+dt*M3,p+dt*N3,activeindex)

        c1=c1+(dt/6)*(K1+2*K2+2*K3+K4)
        c2=c2+(dt/6)*(L1+2*L2+2*L3+L4)
        x=x+(dt/6)*(M1+2*M2+2*M3+M4)
        p=p+(dt/6)*(N1+2*N2+2*N3+N4)    

        return
    endsubroutine

!********************************************************
!                  velocity adjustment                  *
!********************************************************
    subroutine v_adjust(x,p,activeindex)
        implicit none
        double precision :: x,p
        integer :: activeindex

        if (activeindex==1) then
            p=sqrt(1+(2*m*(E1(x)-E2(x)))/(p**2))*p
        elseif (activeindex==2) then
            p=sqrt(1+(2*m*(E2(x)-E1(x)))/(p**2))*p
        endif

        return
    endsubroutine

!********************************************************
!                  Branching corrected                  *
!********************************************************
    subroutine BranchCorrect(t,x,p,c1,c2,activeindex)
        implicit none
        double precision :: t,x,p
        complex(kind=8) :: c1,c2
        integer :: activeindex

        double precision :: pna

        if (activeindex==1) then
            if ((p**2/(2*m)+E1(x)-E2(x))>0) then
                pna=sqrt(1+(2*m*(E1(x)-E2(x)))/(p**2))*p
            elseif ((p**2/(2*m)+E1(x)-E2(x))<0) then
                pna=(p/abs(p))*1.0D-5
            endif
            
            if ((p*(-D_E1(x))*(p-D_E1(x)*dt)*(-D_E1(x)))<0) then
                c2=0
                c1=c1/abs(c1)
            elseif ((pna*(-D_E2(x))*(pna-D_E2(x)*dt)*(-D_E2(x)))<0) then
                c2=0
                c1=c1/abs(c1)
            endif

        elseif (activeindex==2) then
            pna=sqrt(1+(2*m*(E2(x)-E1(x)))/(p**2))*p

            if ((p*(-D_E2(x))*(p-D_E2(x)*dt)*(-D_E2(x)))<0) then
                c2=c2/abs(c2)
                c1=0
            elseif ((pna*(-D_E1(x))*(pna-D_E1(x)*dt)*(-D_E1(x)))<0) then
                c2=c2/abs(c2)
                c1=0
            endif
        endif

        return 
    endsubroutine

!********************************************************
!                 evolution of the system               *
!********************************************************
    subroutine integrationofmotion(t,c1,c2,x,p,rho,activeindex)
        implicit none
        double precision :: t,x,p
        complex(kind=8) :: c1
        complex(kind=8) :: c2
        complex(kind=8) :: rho(2,2)
        integer :: activeindex

        double precision :: rho0_1,rho0_2
        double precision :: d_rho0_1,d_rho0_2
        double precision :: zeta

        rho0_1=dble(c1*conjg(c1))
        rho0_2=dble(c2*conjg(c2))
        d_rho0_1=-2*real(conjg(c1)*c2*(p/m)*d12(x))
        d_rho0_2=-2*real(conjg(c2)*c1*(p/m)*d21(x))

        call RK4(t,c1,c2,x,p,activeindex)

        !call random_seed
        call random_number(zeta)

        if (activeindex==1) then
            if ((zeta<(-dt*d_rho0_1/rho0_1)) .and. (d_rho0_1<0).and. (p**2/(2*m)+E1(x)-E2(x)>0)) then
                call v_adjust(x,p,activeindex)
                activeindex=2
            endif
        elseif (activeindex==2) then
            if (zeta<(-dt*d_rho0_2/rho0_2) .and. (d_rho0_2<0)) then
                call v_adjust(x,p,activeindex)
                activeindex=1
            endif
        endif

        call BranchCorrect(t,x,p,c1,c2,activeindex)

        rho(1,1)=c1*conjg(c1)
        rho(1,2)=c1*conjg(c2)
        rho(2,1)=c2*conjg(c1)
        rho(2,2)=c2*conjg(c2)

        return 
    endsubroutine

end module dynamics


program main
    use constant
    use Hamiltonian
    use dynamics
    implicit none
    double precision :: x
    double precision :: p
    double precision :: x0
    double precision :: p0
    double precision :: t=0
    double precision :: xleft,xfinal
    double precision :: g1,g2
    double precision :: sigma_x
    double precision :: sigma_p
    integer :: Ntrj=1000
    double precision :: dx=0.1
    integer :: activeindex
    complex(kind=8) :: c1,c2
    complex(kind=8) :: rho(2,2)

    double precision :: E_total
    integer :: itrj=1
    double precision :: Reflectlower,Transmitupper,Transmitlower,Reflectupper
    
    Reflectlower=0
    Transmitupper=0
    Transmitlower=0
    Reflectupper=0

    !open(1,file='adiabaticPES.txt')
    !open(2,file='population.txt')
    !open(3,file='nuclearenergy.txt')
    !open(4,file='scattering.txt')

    !write(4,*) 'Reflectlower    ','Reflectupper    ','Transmitlower    ','Transmitupper'

    ! do while(x<10)
    !     write(1,*)x,E1(x),E2(x),d12(x)/50,D_E1(x),D_E2(x),(E1(x+dx)-E1(x))/dx,(E2(x+dx)-E2(x))/dx
    !     x=x+dx
    ! enddo

    read(*,*) p0
    xleft=-17.5
    sigma_x=10/p0
    sigma_p=1/(2*sigma_x)
    x0=xleft-3*sigma_x
    xfinal=33.0

    call random_seed()

    do itrj=1,Ntrj
        call Gaussian_generate(g1,g2)
        x=x0+sigma_x*g1
        p=p0+sigma_p*g2
        c1=1
        c2=0
        activeindex=1

        rho(1,1)=c1*conjg(c1)
        rho(1,2)=c1*conjg(c2)
        rho(2,1)=c2*conjg(c1)
        rho(2,2)=c2*conjg(c2)

        do while(x<=xfinal .and. x>=-xfinal)
            call integrationofmotion(t,c1,c2,x,p,rho,activeindex)
            t=t+dt
            ! if (activeindex==1) then
            !     E_total=E1(x)
            ! elseif (activeindex==2) then
            !     E_total=E2(x)
            ! endif

            ! write(2,*)t,dble(c1*conjg(c1)),dble(c2*conjg(c2)),dble(c1*conjg(c1)+c2*conjg(c2))
            ! write(3,*)t,p**2/(2*m),E_total,p**2/(2*m)+E_total
        enddo

        if ((x>xfinal) .and. (activeindex==1)) then
            Transmitlower=Transmitlower+1
        elseif ((x>xfinal) .and. (activeindex==2)) then
            Transmitupper=Transmitupper+1
        elseif ((x<-xfinal) .and. (activeindex==1)) then
            Reflectlower=Reflectlower+1
        elseif ((x<-xfinal) .and. (activeindex==2)) then
            Reflectupper=Reflectupper+1
        endif
    enddo

    Reflectlower=Reflectlower/Ntrj
    Transmitupper=Transmitupper/Ntrj
    Transmitlower=Transmitlower/Ntrj
    Reflectupper=Reflectupper/Ntrj

    write(*,"(F6.3,4X,F12.9,4X,F12.9,4X,F12.9,4X,F12.9)") p0,Reflectlower,Reflectupper,Transmitlower,Transmitupper

    !close(1)
    !close(2)
    !close(3)
    !close(4)
    stop
end program main

subroutine Gaussian_generate(g1,g2)
    use constant
    implicit none
    double precision :: g1,g2
    double precision :: x1,x2

    call random_number(x1)
    call random_number(x2)
    
    g1=sqrt(-2*log(x1))*cos(2*pi*x2)
    g2=sqrt(-2*log(x1))*sin(2*pi*x2)

    return
endsubroutine
