!-----------------------------------------------------------------------
!     PMLの共通変数
!-----------------------------------------------------------------------
      module fdtd_pml
      use fdtd_variable 

      type pml_var
         integer::nx0,ny0,nz0,nx1,ny1,nz1
         real,pointer::exy(:,:,:),exz(:,:,:)
         real,pointer::eyz(:,:,:),eyx(:,:,:)
         real,pointer::ezx(:,:,:),ezy(:,:,:)
         real,pointer::hxy(:,:,:),hxz(:,:,:)
         real,pointer::hyz(:,:,:),hyx(:,:,:)
         real,pointer::hzx(:,:,:),hzy(:,:,:)
         real,pointer::aeyx(:,:,:),amyx(:,:,:)
         real,pointer::aezx(:,:,:),amzx(:,:,:)
         real,pointer::aexy(:,:,:),amxy(:,:,:)
         real,pointer::aezy(:,:,:),amzy(:,:,:)
         real,pointer::aexz(:,:,:),amxz(:,:,:)
         real,pointer::aeyz(:,:,:),amyz(:,:,:)
         real,pointer::beyx(:,:,:),bmyx(:,:,:)
         real,pointer::bezx(:,:,:),bmzx(:,:,:)
         real,pointer::bexy(:,:,:),bmxy(:,:,:)
         real,pointer::bezy(:,:,:),bmzy(:,:,:)
         real,pointer::bexz(:,:,:),bmxz(:,:,:)
         real,pointer::beyz(:,:,:),bmyz(:,:,:)
      end type pml_var

      type(pml_var)::pml_x0,pml_x1,pml_y0,pml_y1,pml_z0,pml_z1

      real,parameter::copml=-1.5280063e-4
      end module fdtd_pml
!-----------------------------------------------------------------------
!     初期設定
!-----------------------------------------------------------------------
      subroutine init_pml()
      use fdtd_pml

      call addpml(pml_x0,       0, lpml,       0,   ny,       0,   nz) 
      call addpml(pml_x1, nx-lpml,   nx,       0,   ny,       0,   nz) 
      call addpml(pml_y0,       0,   nx,       0, lpml,       0,   nz) 
      call addpml(pml_y1,       0,   nx, ny-lpml,   ny,       0,   nz) 
      call addpml(pml_z0,       0,   nx,       0,   ny,       0, lpml) 
      call addpml(pml_z1,       0,   nx,       0,   ny, nz-lpml,   nz)
      end subroutine init_pml
!-----------------------------------------------------------------------
!     電界に対するPML
!-----------------------------------------------------------------------
      subroutine e_pml()
      use fdtd_pml

      call epml(pml_x0)
      call epml(pml_x1) 
      call epml(pml_y0)
      call epml(pml_y1) 
      call epml(pml_z0)
      call epml(pml_z1) 
      end subroutine e_pml
!-----------------------------------------------------------------------
!     磁界に対するPML
!-----------------------------------------------------------------------
      subroutine h_pml()
      use fdtd_pml

      call hpml(pml_x0)
      call hpml(pml_x1) 
      call hpml(pml_y0)
      call hpml(pml_y1) 
      call hpml(pml_z0)
      call hpml(pml_z1)
      end subroutine h_pml
!-----------------------------------------------------------------------
!     係数の計算
!-----------------------------------------------------------------------
      subroutine addpml(pml,nx0,nx1,ny0,ny1,nz0,nz1)
      use fdtd_pml
      type(pml_var)::pml
      real::mupml

      pml%nx0=nx0
      pml%nx1=nx1
      pml%ny0=ny0
      pml%ny1=ny1
      pml%nz0=nz0
      pml%nz1=nz1
      allocate(pml%exy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%exz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%eyz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%eyx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%ezx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%ezy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%hxy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%hxz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%hyz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%hyx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%hzx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%hzy(nx0:nx1,ny0:ny1,nz0:nz1))

      allocate(pml%aeyx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%aezx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%aexy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%aezy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%aexz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%aeyz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%amyx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%amzx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%amxy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%amzy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%amxz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%amyz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%beyx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bezx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bexy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bezy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bexz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%beyz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bmyx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bmzx(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bmxy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bmzy(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bmxz(nx0:nx1,ny0:ny1,nz0:nz1))
      allocate(pml%bmyz(nx0:nx1,ny0:ny1,nz0:nz1))
      pml%exy=0.0
      pml%exz=0.0
      pml%eyz=0.0
      pml%eyx=0.0
      pml%ezx=0.0
      pml%ezy=0.0
      pml%hxy=0.0
      pml%hxz=0.0
      pml%hyz=0.0
      pml%hyx=0.0
      pml%hzx=0.0
      pml%hzy=0.0

      smax0x=copml*rmax*(order+1)/(lpml*dx)
      smax0y=copml*rmax*(order+1)/(lpml*dy)
      smax0z=copml*rmax*(order+1)/(lpml*dz)
 
      do k=nz0,nz1-1
         do j=ny0,ny1-1
            do i=nx0,nx1-1
             if(i<lpml) then
                sigmxm=((lpml-i-0.5)/lpml)**order*smax0x
                sigmxe=(float(lpml-i)/lpml)**order*smax0x
             else if(i>=nx-lpml) then
                sigmxm=((i-nx+lpml+0.5)/lpml)**order*smax0x
                sigmxe=(float(i-nx+lpml)/lpml)**order*smax0x
             else                                         
                sigmxm=0.0
                sigmxe=0.0
             end if

             if(j<lpml) then
                sigmym=((lpml-j-0.5)/lpml)**order*smax0y
                sigmye=(float(lpml-j)/lpml)**order*smax0y
             else if(j>=ny-lpml) then
                sigmym=((j-ny+lpml+0.5)/lpml)**order*smax0y
                sigmye=(float(j-ny+lpml)/lpml)**order*smax0y
             else                                         
                sigmym=0.0
                sigmye=0.0
             end if
           
             if(k<lpml) then
                sigmzm=((lpml-k-0.5)/lpml)**order*smax0z
                sigmze=(float(lpml-k)/lpml)**order*smax0z
             else if(k>=nz-lpml) then
                sigmzm=((k-nz+lpml+0.5)/lpml)**order*smax0z
                sigmze=(float(k-nz+lpml)/lpml)**order*smax0z
             else                                         
                sigmzm=0.0
                sigmze=0.0
             end if

!Exの係数
             epspml=0.25*(epsd(i,j,k)+epsd(i,j-1,k)
     &              +epsd(i,j,k-1)+epsd(i,j-1,k-1))*eps0
             sigmy=sigmye*(epspml/eps0)
             sigmz=sigmze*(epspml/eps0)

             a=0.5*sigmy*dt/epspml
             pml%aexy(i,j,k)=(1.0-a)/(1.0+a)
             pml%bexy(i,j,k)=dt/epspml/(1.0+a)/dy
c
             a=0.5*sigmz*dt/epspml
             pml%aexz(i,j,k)=(1.0-a)/(1.0+a)
             pml%bexz(i,j,k)=dt/epspml/(1.0+a)/dz

!Eyの係数
             epspml=0.25*(epsd(i,j,k)+epsd(i-1,j,k)
     &              +epsd(i,j,k-1)+epsd(i-1,j,k-1))*eps0
             sigmz=sigmze*(epspml/eps0)
             sigmx=sigmxe*(epspml/eps0)

             a=0.5*sigmz*dt/epspml
             pml%aeyz(i,j,k)=(1.0-a)/(1.0+a)
             pml%beyz(i,j,k)=dt/epspml/(1.0+a)/dz
c
             a=0.5*sigmx*dt/epspml
             pml%aeyx(i,j,k)=(1.0-a)/(1.0+a)
             pml%beyx(i,j,k)=dt/epspml/(1.0+a)/dx

!Ezの係数
             epspml=0.25*(epsd(i,j,k)+epsd(i-1,j,k)
     &              +epsd(i,j-1,k)+epsd(i-1,j-1,k))*eps0
             sigmx=sigmxe*(epspml/eps0)
             sigmy=sigmye*(epspml/eps0)

             a=0.5*sigmx*dt/epspml
             pml%aezx(i,j,k)=(1.0-a)/(1.0+a)
             pml%bezx(i,j,k)=dt/epspml/(1.0+a)/dx
c
             a=0.5*sigmy*dt/epspml
             pml%aezy(i,j,k)=(1.0-a)/(1.0+a)
             pml%bezy(i,j,k)=dt/epspml/(1.0+a)/dy

!Hxの係数
             mupml=0.5*(mud(i,j,k)+mud(i-1,j,k))*mu0
             epspml=0.5*(epsd(i,j,k)+epsd(i-1,j,k))*eps0
             sigmy=sigmym*(epspml/eps0)
             sigmz=sigmzm*(epspml/eps0)

             a=0.5*sigmy*dt/epspml
             pml%amxy(i,j,k)=(1.0-a)/(1.0+a)
             pml%bmxy(i,j,k)=dt/mupml/(1.0+a)/dy
c
             a=0.5*sigmz*dt/epspml
             pml%amxz(i,j,k)=(1.0-a)/(1.0+a)
             pml%bmxz(i,j,k)=dt/mupml/(1.0+a)/dz

!Hyの係数
             mupml=0.5*(mud(i,j,k)+mud(i,j-1,k))*mu0
             epspml=0.5*(epsd(i,j,k)+epsd(i,j-1,k))*eps0
             sigmz=sigmzm*(epspml/eps0)
             sigmx=sigmxm*(epspml/eps0)

             a=0.5*sigmz*dt/epspml
             pml%amyz(i,j,k)=(1.0-a)/(1.0+a)
             pml%bmyz(i,j,k)=dt/mupml/(1.0+a)/dz
c
             a=0.5*sigmx*dt/epspml
             pml%amyx(i,j,k)=(1.0-a)/(1.0+a)
             pml%bmyx(i,j,k)=dt/mupml/(1.0+a)/dx

!Hzの係数
             mupml=0.5*(mud(i,j,k)+mud(i,j,k-1))*mu0
             epspml=0.5*(epsd(i,j,k)+epsd(i,j,k-1))*eps0
             sigmx=sigmxm*(epspml/eps0)
             sigmy=sigmym*(epspml/eps0)

             a=0.5*sigmx*dt/epspml
             pml%amzx(i,j,k)=(1.0-a)/(1.0+a)
             pml%bmzx(i,j,k)=dt/mupml/(1.0+a)/dx
c
             a=0.5*sigmy*dt/epspml
             pml%amzy(i,j,k)=(1.0-a)/(1.0+a)
             pml%bmzy(i,j,k)=dt/mupml/(1.0+a)/dy
           end do
         end do
      end do               
      end subroutine addpml
!-----------------------------------------------------------------------
!     電界の計算
!-----------------------------------------------------------------------
      subroutine epml(pml)
      use fdtd_pml
      type(pml_var)::pml
!Ex
      do k=pml%nz0+1,pml%nz1-1
         do j=pml%ny0+1,pml%ny1-1
            do i=pml%nx0,pml%nx1-1
               pml%exy(i,j,k)=pml%aexy(i,j,k)*pml%exy(i,j,k)
     &               +pml%bexy(i,j,k)*(hz(i,j,k)-hz(i,j-1,k))
               pml%exz(i,j,k)=pml%aexz(i,j,k)*pml%exz(i,j,k)
     &               +pml%bexz(i,j,k)*(hy(i,j,k-1)-hy(i,j,k))
               ex(i,j,k)=pml%exy(i,j,k)+pml%exz(i,j,k)
            end do
         end do
      end do
!Ey
      do k=pml%nz0+1,pml%nz1-1
         do j=pml%ny0,pml%ny1-1
            do i=pml%nx0+1,pml%nx1-1
               pml%eyz(i,j,k)=pml%aeyz(i,j,k)*pml%eyz(i,j,k)
     &               +pml%beyz(i,j,k)*(hx(i,j,k)-hx(i,j,k-1))
               pml%eyx(i,j,k)=pml%aeyx(i,j,k)*pml%eyx(i,j,k)
     &               +pml%beyx(i,j,k)*(hz(i-1,j,k)-hz(i,j,k))
               ey(i,j,k)=pml%eyz(i,j,k)+pml%eyx(i,j,k)
           end do
         end do
      end do
!Ez
      do k=pml%nz0,pml%nz1-1
         do j=pml%ny0+1,pml%ny1-1
            do i=pml%nx0+1,pml%nx1-1
               pml%ezx(i,j,k)=pml%aezx(i,j,k)*pml%ezx(i,j,k)
     &               +pml%bezx(i,j,k)*(hy(i,j,k)-hy(i-1,j,k))
               pml%ezy(i,j,k)=pml%aezy(i,j,k)*pml%ezy(i,j,k)
     &               +pml%bezy(i,j,k)*(hx(i,j-1,k)-hx(i,j,k))
               ez(i,j,k)=pml%ezx(i,j,k)+pml%ezy(i,j,k)
            end do
         end do
      end do
      end subroutine epml
!-----------------------------------------------------------------------
!     磁界の計算
!-----------------------------------------------------------------------
      subroutine hpml(pml)      
      use fdtd_pml
      type(pml_var)::pml
!Hx
      do k=pml%nz0,pml%nz1-1
         do j=pml%ny0,pml%ny1-1
            do i=pml%nx0+1,pml%nx1-1
               pml%hxy(i,j,k)=pml%amxy(i,j,k)*pml%hxy(i,j,k)
     &               +pml%bmxy(i,j,k)*(ez(i,j,k)-ez(i,j+1,k))
               pml%hxz(i,j,k)=pml%amxz(i,j,k)*pml%hxz(i,j,k)
     &               +pml%bmxz(i,j,k)*(ey(i,j,k+1)-ey(i,j,k))
               hx(i,j,k)=pml%hxy(i,j,k)+pml%hxz(i,j,k)
            end do
         end do
      end do
!Hy
      do k=pml%nz0,pml%nz1-1
         do j=pml%ny0+1,pml%ny1-1
            do i=pml%nx0,pml%nx1-1
               pml%hyz(i,j,k)=pml%amyz(i,j,k)*pml%hyz(i,j,k)
     &               +pml%bmyz(i,j,k)*(ex(i,j,k)-ex(i,j,k+1))
               pml%hyx(i,j,k)=pml%amyx(i,j,k)*pml%hyx(i,j,k)
     &               +pml%bmyx(i,j,k)*(ez(i+1,j,k)-ez(i,j,k))
               hy(i,j,k)=pml%hyz(i,j,k)+pml%hyx(i,j,k)
            end do
         end do
      end do
!Hz
      do k=pml%nz0+1,pml%nz1-1
         do j=pml%ny0,pml%ny1-1
            do i=pml%nx0,pml%nx1-1
               pml%hzx(i,j,k)=pml%amzx(i,j,k)*pml%hzx(i,j,k)
     &               +pml%bmzx(i,j,k)*(ey(i,j,k)-ey(i+1,j,k))
               pml%hzy(i,j,k)=pml%amzy(i,j,k)*pml%hzy(i,j,k)
     &               +pml%bmzy(i,j,k)*(ex(i,j+1,k)-ex(i,j,k))
               hz(i,j,k)=pml%hzx(i,j,k)+pml%hzy(i,j,k)
            end do
         end do
      end do 
      end subroutine hpml
