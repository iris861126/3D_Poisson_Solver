! use uft-8 in Linux, use ANSI in Windows
! compile
! f2py -m poisson_ter -c poisson_ter.f95 --opt=-O3

function check_bdy_type(bdy_type)
implicit none
! 新增地形的處理方式
! 遇到地形 那一個就是下邊界了
! 再下面就跳過不疊代
character*8 :: bdy_type
integer :: check_bdy_type

if (trim(bdy_type) == 'zero') then
    check_bdy_type = 0
elseif (trim(bdy_type) == 'Neumann') then
    check_bdy_type = 1
else
    write(*,*) 'Unknown Boundary type ' // trim(bdy_type) // '.'
    stop 1
endif

return
end function

subroutine poisson3D_solver_SOR(INI,RHS,BCx,BCy,BCz,HGT_index,dx,dy,dz,iter_max,tol,omega,x_bdy_type,y_bdy_type,z_bdy_type &
                              	&,Nx,Ny,Nz,psiout,errorout,iter_count_out)
implicit none

integer :: check_bdy_type

!intend in
real*8, intent(in) :: dx, dy, dz
integer, intent(in) :: Nx, Ny, Nz
real*8, intent(in) :: BCx(Nz,Ny,Nx), BCy(Nz,Ny,Nx), BCz(Nz,Ny,Nx), RHS(Nz,Ny,Nx), HGT_index(Ny,Nx), INI(2,Nz,Ny,Nx)
integer, intent(in) :: iter_max
real*8, intent(in) :: tol, omega
character*8, intent(in) :: x_bdy_type, y_bdy_type, z_bdy_type

!intent out
real*8 :: psi(2,Nz,Ny,Nx)
real*8, intent(out) :: psiout(Nz,Ny,Nx)
real*8, intent(out) :: errorout
integer, intent(out) :: iter_count_out

!state index
integer :: new = 2-1, old = 1-1
integer :: t, i, j, k     !iteration idx
integer :: x_bdy, y_bdy, z_bdy

!others
real*8 :: h, a, b, c, d, error = 1.

! 事前準備
psi = INI

x_bdy = check_bdy_type(x_bdy_type)
y_bdy = check_bdy_type(y_bdy_type)
z_bdy = check_bdy_type(z_bdy_type)

h = dy*dy*dz*dz + dx*dx*dz*dz + dx*dx*dy*dy
a = dy*dy*dz*dz/2/h
b = dx*dx*dz*dz/2/h
c = dx*dx*dy*dy/2/h
d = dx*dx*dy*dy*dz*dz/2/h


! 迭代
do t = 1, iter_max
    old = mod(old,2)+1
    new = mod(new,2)+1
    
    do i = 1, Nx
        do j = 1, Ny
            do k = 1, Nz
                !處理邊界
                if (k < HGT_index(j,i)+1) then   !在地形以下
                    cycle                      	 !跳過不處理
                elseif (i == 1) then !-x
                    if (x_bdy == 0) then
                        psi(new,k,j,i) = 0
                    elseif (x_bdy == 1) then
                        psi(new,k,j,i) = psi(old,k,j,i+2) - 2*dx*BCx(k,j,i+1)
                    endif
                elseif (i == Nx) then !+x
                    if (x_bdy == 0) then
                        psi(new,k,j,i) = 0
                    elseif (x_bdy == 1) then
                        psi(new,k,j,i) = psi(old,k,j,i-2) + 2*dx*BCx(k,j,i-1)
                    endif
                elseif (j == 1) then !-y
                    if (y_bdy == 0) then
                        psi(new,k,j,i) = 0
                    elseif (y_bdy == 1) then
                        psi(new,k,j,i) = psi(old,k,j+2,i) - 2*dy*BCy(k,j+1,i)
                    endif
                elseif (j == Ny) then !+y
                    if (y_bdy == 0) then
                        psi(new,k,j,i) = 0
                    elseif (y_bdy == 1) then
                        psi(new,k,j,i) = psi(old,k,j-2,i) + 2*dy*BCy(k,j-1,i)
                    endif
                elseif ((k == 1) .or. (k == HGT_index(j,i)+1)) then !-z 在網格底部、在最靠近地形的點都是下邊界
                    psi(new,k,j,i) = psi(old,k+2,j,i) - 2*dz*BCz(k+1,j,i)
                elseif (k == Nz) then !+z
                    if (z_bdy == 0) then
                        psi(new,k,j,i) = 0
                    elseif (z_bdy == 1) then
                        psi(new,k,j,i) = psi(old,k-2,j,i) + 2*dz*BCz(k-1,j,i)
                    endif
                else
             
                    !內部點計算
                    psi(new,k,j,i) = (a*(psi(new,k,j,i-1) + psi(old,k,j,i+1)) &
                                 & + b*(psi(new,k,j-1,i) + psi(old,k,j+1,i)) &
                                 & + c*(psi(new,k-1,j,i) + psi(old,k+1,j,i)) &
                                 & - d*RHS(k,j,i))*omega + (1-omega)*psi(old,k,j,i)
                endif
            end do
        end do
    end do
    !write(*,*) psi(new,3,5,5)
    !收斂判定
    if (mod(t,10) == 1) then
        error = MAXVAL(dabs(psi(1,:,:,:)-psi(2,:,:,:)))/omega
        write(*,*) t,error
    end if  
    if (error < tol) then
        exit
    end if
end do
error = MAXVAL(dabs(psi(1,:,:,:)-psi(2,:,:,:)))/omega

!邊點處理
do i = 1, Nx 
    psi(new,1,1,i) = (psi(new,2,1,i) + psi(new,1,2,i))/2
    psi(new,Nz,1,i) = (psi(new,Nz-1,1,i) + psi(new,Nz,2,i))/2
    psi(new,1,Ny,i) = (psi(new,2,Ny,i) + psi(new,1,Ny-1,i))/2
    psi(new,Nz,Ny,i) = (psi(new,Nz-1,Ny,i) + psi(new,Nz,Ny-1,i))/2
end do

do j = 1, Ny
    psi(new,1,j,1) = (psi(new,2,j,1) + psi(new,1,j,2))/2
    psi(new,1,j,Nx) = (psi(new,2,j,Nx) + psi(new,1,j,Nx-1))/2
    psi(new,Nz,j,1) = (psi(new,Nz-1,j,1) + psi(new,Nz,j,2))/2
    psi(new,Nz,j,Nx) = (psi(new,Nz-1,j,Nx) + psi(new,Nz,j,Nx-1))/2
end do

do k = 1, Nz
    psi(new,k,1,1) = (psi(new,k,2,1) + psi(new,k,1,2))/2
    psi(new,k,Ny,1) = (psi(new,k,Ny-1,1) + psi(new,k,Ny,2))/2
    psi(new,k,1,Nx) = (psi(new,k,2,Nx) + psi(new,j,k,Nx-1))/2
    psi(new,k,Ny,Nx) = (psi(new,k,Ny-1,Nx) + psi(new,k,Ny,Nx-1))/2
end do

!角落點處理
psi(new,1,1,1)    = (psi(new,2,1,1) + psi(new,1,2,1) + psi(new,1,1,2))/3
psi(new,Nz,1,1)   = (psi(new,Nz-1,1,1) + psi(new,Nz,2,1) + psi(new,Nz,1,2))/3
psi(new,1,Ny,1)   = (psi(new,2,Ny,1) + psi(new,1,Ny-1,1) + psi(new,1,Ny,2))/3
psi(new,1,1,Nx)   = (psi(new,2,1,Nx) + psi(new,1,2,Nx) + psi(new,1,1,Nx-1))/3
psi(new,Nz,Ny,1)  = (psi(new,Nz-1,Ny,1) + psi(new,Nz,Ny-1,1) + psi(new,Nz,Ny,2))/3
psi(new,1,Ny,Nx)  = (psi(new,2,Ny,Nx) + psi(new,1,Ny-1,Nx) + psi(new,1,Ny,Nx-1))/3
psi(new,Nz,1,Nx)  = (psi(new,Nz-1,1,Nx) + psi(new,Nz,2,Nx) + psi(new,Nz,1,Nx-1))/3
psi(new,Nz,Ny,Nx) = (psi(new,Nz-1,Ny,Nx) + psi(new,Nz,Ny-1,Nx) + psi(new,Nz,Ny,Nx-1))/3

!輸出處理
psiout = psi(new,:,:,:)
errorout = error
iter_count_out = t

end subroutine