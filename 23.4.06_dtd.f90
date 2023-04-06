	!********************************START******************************
	program main
	implicit none
	!********************************����*******************************
	!�߽���ڲ�����
	real(kind=8):: xs(2,2,201),ys(2,2,201)
	real(kind=8),allocatable :: zlast(:,:,:,:),xlast(:,:,:,:),ylast(:,:,:,:)
	real(kind=8),allocatable :: zlast1(:,:,:,:),xlast1(:,:,:,:),ylast1(:,:,:,:)
	!�ֲ�����
	real(kind=8) :: s(301),sb(2,301)
	!�ߴ�
	real(kind=8):: theta,theta3,theta4
	real(kind=8),parameter:: pi=4.0*atan(1.0)
	real(kind=8) :: zd2,temp1,temp3
	!���������Կ��ƺ�����ѭ������
	real(kind=8) :: p,q
	integer :: n,j,m,k,i
	!Z����ľ���
	!real(kind=8) :: lz
	!real(kind=8) :: zjl(301,301)
	!-----------------��������������--------------------------------------
	!�ڲ�����
	real(kind=8):: yt(301,301),zt(301,301),xt(301,301)
	! �߽�����
	real(kind=8):: xss(2,301),yss(2,301),zss(2,301)
	real(kind=8):: ymm(301,301,25),zmm(301,301,25)
	!�ڵ���
	integer :: nzt(31),nxt(31),nyt(31),imax(31),jmax(31),kmax(31)
	!PLOT-3D��ʽ�ڵ���������ɫ����
	integer::  nblocks,mf
	!����
	real(kind=8):: tempxp,tempyp
	!�ڵ���
	integer:: nytemp,nxtemp,nxra,nyrb,jj1,jj2,jj3,jj4
	integer:: nra,nrb,na,mb,nii,mjj
	!�ߴ�
	real(kind=8):: ra,rb,rh,alfa,byta
	!��ʱһά����
	real(kind=8):: temporary(20)
	!��ʱ��ά����
	real(kind=8):: temporary1(15,301),temporary2(8,301),temporary3(6,35)
	!------------ǿ��������--------------------------------------------------------
	real(kind=8):: ztnew1(301,301),xtnew1(301,301),ytnew1(301,301)
	real(kind=8):: rt(301,201),thetat(301,201)
	real(kind=8):: de,xL,da,yL
	!------------����Բ��---------------------------
	real(kind=8):: ztnew2(301,301),xtnew2(301,301),ytnew2(301,301)
	real(kind=8):: dm
	!********************************�ߴ�*******************************
	xL=4.47046286
	yL=1.10993419
	de=xL/(pi)
	da=0.6312*de
	nra=8
	nrb=3
	ra=xl/float(nra)
	rb=yl/float(nrb)
	alfa=atan(rb/ra)
	!���²�����������
	rh=0.03557481 !��������ĸ߶�
	dm=1.51 !��Բ��ֱ��
	zd2=0.5
	temp1=0.5
	temp3=0.5
	!lz=30.0
	!********************************����������***************************
	!�������ֲ�Ԫģ�ͽ���-˫�߽編
	!��X������м��ܴ���
	nxtemp=9
	n=int(nxtemp/2)+1
	p=0.64585
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	!д���±߽���������
	do i=1,2*n-1
		xss(1,i)=(ra/2.0)*sb(1,i)
		yss(1,i)=0.0
		xss(2,i)=xss(1,i)
		yss(2,i)=rb/2.0
	end do
	!--------------------------------------------------------------------
	!��Y������м��ܴ���
	nytemp=11
	n=int(nytemp/2)+1
	p=0.64585
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do   j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	!ʹ��˫�߽編д�ڲ�����
	do i=1,nxtemp
		do j=1,nytemp
			xt(i,j)=xss(1,i)*(1.0-sb(2,j))+sb(2,j)*xss(2,i)
			yt(i,j)=yss(1,i)*(1.0-sb(2,j))+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	!����ֲ�Ԫ����
	open(2,file='������-�ֲ�Ԫ����-��ʼ.plt')
	WRITE(2,*) 'variables= "x","y"'
	WRITE(2,*) 'zone i=',nytemp,'j=',nxtemp
	do  i=1,nxtemp
		do   j=1,nytemp
			write(2,*) xt(i,j),yt(i,j)
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!�Ѵֲ�Ԫ��x�����ӳ����������ܹ�ϵ���ֲ��䣩,�ӳ�����Ϊra/2.0
	do  i=1,nxtemp
		do   j=1,nytemp
			xt(nxtemp+i-1,j)=-xt(nxtemp-i+1,j)+ra
			yt(nxtemp+i-1,j)=yt(nxtemp-i+1,j)
		end do
	end do
	!--------------------------------------------------------------------
	!�����x�����ӳ���Ĵֲ�Ԫ����
	open(2,file='������-�ֲ�Ԫ����+X.plt')
	write(2,*) 'variables= "x","y"'
	write(2,*) 'zone i=',nytemp,'j=',nxtemp
	do  i=1,nxtemp
		do   j=1,nytemp
			write(2,*) xt(i,j),yt(i,j)  !ԭ�ȵĲ���
		end do
	end do

	write(2,*) 'zone i=',nytemp,'j=',nxtemp
	do  i=1,nxtemp
		do   j=1,nytemp
			write(2,*) xt(nxtemp+i-1,j),yt(nxtemp+i-1,j)  !�ӳ��Ĳ���
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!�Ѵֲ�Ԫ��Y�����ӳ����������ܹ�ϵ���ֲ��䣩,�ӳ�����Ϊrb/2.0
	nxra=2*nxtemp-1
	do  i=1,nxra
		do   j=1,nytemp
			xt(i,nytemp+j-1)=xt(i,nytemp-j+1)
			yt(i,nytemp+j-1)=-yt(i,nytemp-j+1)+rb
		end do
	end do
	!--------------------------------------------------------------------
	!�����y�����ӳ���Ĵֲ�Ԫ����
	open(2,file='������-�ֲ�Ԫ����+Y.plt')
	WRITE(2,*) 'variables= "x","y"'
	WRITE(2,*) 'zone i=',nytemp,'j=',nxra
	do  i=1,nxra
		do   j=1,nytemp
			write(2,*) xt(i,j),yt(i,j)    !ԭ�ȵĲ���
		end do
	end do
	!--------------------------------------------------------------------
	WRITE(2,*) 'zone i=',nytemp,'j=',nxra
	do  i=1,nxra
		do   j=1,nytemp
			write(2,*) xt(i,nytemp+j-1),yt(i,nytemp+j-1) !�ӳ��Ĳ���
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!�ӳ�֮�������ֲ�Ԫģ��
	nyrb=2*nytemp-1
	open(2,file='������-�ֲ�Ԫ+����.plt')
	WRITE(2,*) 'variables= "x","y"'
	WRITE(2,*) 'zone i=',nyrb,'j=',nxra
	do  i=1,nxra
		do   j=1,nyrb
			write(2,*) xt(i,j),yt(i,j)
		end do
	end do
	close(2)
	!write(*,*)'yt',yt!,���Ϊʲô���ȫ��0��
	!pause
	write(*,*)'�ֲ�Ԫ-����-xt(k,1)',(xt(k,1),k=1,nxra)      !=ra,�ֲ�Ԫ�ĳ�
	pause
	write(*,*)'�ֲ�Ԫ-����-yt(1,k)',(yt(1,k),k=1,nyrb)     !=rb���ֲ�Ԫ�Ŀ�
	pause
	!--------------------------------------------------------------------
	!���潨����������ģ��
	do i=1,nxra
		do j=1,nyrb
			tempxp=abs(xt(i,j))  !���ز����ľ���ֵ
			tempyp=abs(yt(i,j))
			!�����Ǵֲ�Ԫ������������ĳ��򣬵��ǿ�����
			if(tempxp.ge.0.0.and.tempxp.le.0.50*ra.and. &
				&   tempyp.ge.0.0.and.tempyp.le.tempxp*tan(alfa*1.0)) then
				zt(i,j)=2.0*rh/ra*(tempxp-ra/rb*tempyp)
			end if
			if(tempxp.gt.0.50*ra.and.tempxp.le.ra.and. &
				&   tempyp.ge.0.0d0.and.tempyp.le.(ra-tempxp)*tan(alfa*1.0))  then
				zt(i,j)=2.0*rh/ra*(ra-tempxp-ra/rb*tempyp)
			end if
			if(tempxp.ge.0.0.and.tempxp.le.0.50*ra.and.  &
				&   tempyp.gt.tempxp*tan(alfa*1.0).and.tempyp.le.0.50*rb)  then
				zt(i,j)=2.0*rh/rb*(tempyp-rb/ra*tempxp)
			end if
			if(tempxp.gt.0.50*ra.and.tempxp.le.ra.and.  &
				&   tempyp.gt.(ra-tempxp)*tan(alfa*1.0).and.tempyp.le.0.50*rb)  then
				zt(i,j)=2.0*rh/rb*(tempyp-rb+rb/ra*tempxp)
			end if

			if(tempxp.ge.0.0.and.tempxp.le.0.50*ra.and. &
				&   tempyp.ge.0.50*rb.and.tempyp.le.0.50*rb   &
				&   +(0.50*ra-tempxp)*tan(alfa*1.0))   then
				zt(i,j)=2.0*rh/rb*(-rb/ra*tempxp-tempyp+rb)
			end if
			if(tempxp.ge.0.0.and.tempxp.le.0.50*ra.and. &
				&   tempyp.gt.0.50*rb+(0.50*ra-tempxp)*tan(alfa*1.0).and. &
				&   tempyp.le.rb)   then
				zt(i,j)=2.0*rh/rb*(rb/ra*tempxp+tempyp-rb)
			end if
			if(tempxp.gt.0.50*ra.and.tempxp.le.ra.and. &
				&   tempyp.ge.0.50*rb+(tempxp-0.50*ra)*tan(alfa*1.0).and. &
				&   tempyp.le.rb)    then
				zt(i,j)=2.0*rh/rb*(-rb/ra*tempxp+tempyp)
			end if
			if(tempxp.gt.0.50*ra.and.tempxp.le.ra.and. &
				&   tempyp.gt.0.50*rb.and.tempyp.le.0.50*rb &
				&   +(tempxp-0.50*ra)*tan(alfa*1.0))   then
				zt(i,j)=2.0*rh/ra*(tempxp-ra/rb*tempyp)
			end if
		end do
	end do
	!--------------------------------------------------------------------
	open(2,file='��������-��ʼ.plt')
	WRITE(2,*) 'variables="x", "y","z"'
	WRITE(2,*) 'zone i=',2,'j=',nyrb,'k=',nxra   !�����2����y,z���������ÿ���ظ�2
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,2
				write(2,*) xt(i,j),yt(i,j),zt(i,j)
			end do
		end do
	end do
	close(2)
	write(*,*)'��������-��ʼ-xt(k,1)',(xt(k,1),k=1,nxra)      !=ra,��������ĳ�
	pause
	write(*,*)'��������-��ʼ-yt(1,k)',(yt(1,k),k=1,nyrb)   !=rb,��������Ŀ�
	pause
	!--------------------------------------------------------------------
	!���������������չ
	do  na=1,nra
		do mb=1,nrb
			nii=(na-1)*(nxra-1)
			mjj=(mb-1)*(nyrb-1)
			do i=1,nxra
				do j=1,nyrb
					xt(nii+i,mjj+j)=xt(i,j)+(na-1)*ra
					yt(nii+i,mjj+j)=yt(i,j)+(mb-1)*rb
					zt(nii+i,mjj+j)=zt(i,j)
				end  do
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!�������������������ģ�͸��Ƶ�һ����λ�ã����������Ϊ(nra-1)*ra��(nrb-1)*rb��
	open(2,file='�������������Ͻ���չ.plt')
	WRITE(2,*) 'variables="x", "y","z"'
	WRITE(2,*) 'zone i=',2,'j=',nyrb,'k=',nxra
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,2
				write(2,*)  xt(nii+i,mjj+j),yt(nii+i,mjj+j),zt(nii+i,mjj+j)
			end do            !��չ��Ĳ���
		end do
	end do
	!--------------------------------------------------------------------
	WRITE(2,*) 'zone i=',2,'j=',nyrb,'k=',nxra
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,2
				write(2,*)   xt(i,j),yt(i,j),zt(i,j)
			end do            !��չǰ�Ĳ���
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!��X��Y����Ľڵ������¸�ֵ
	nxtemp=(nra-1)*(nxra-1)+nxra
	nytemp=(nrb-1)*(nyrb-1)+nyrb
	write(6,*)'��չ����½ڵ���nxtemp��nytemp�ֱ��ǣ�' ,nxtemp,nytemp
	pause
	!--------------------------------------------------------------------
	!do i=1,nxtemp
	!	do j=1,nytemp
	!		zjl(i,j)=lz/(nytemp-1)*(j-1)
	!	end do
	!end do

	open(2,file='��չ��Ľ�������.plt')
	WRITE(2,*) 'variables="x", "y","z"'
	WRITE(2,*) 'zone i=',2,'j=',nytemp,'k=',nxtemp
	do i=1,nxtemp			 !��չ��Ľ������峤=nra*ra=xL,��=nrb*rb=yL
		do  j=1,nytemp
			do   k=1,2
				write(2,*)xt(i,j),yt(i,j),zt(i,j)
			end do
		end do
	end do
	close(2)
	write(*,*)'��������-��չ-yt(1,k)',(yt(1,k),k=1,nytemp)
	pause
	!write(*,*)'��չ����������yt(5,k)����', (yt(5,k),k=1,nytemp)
	!write(*,*)'��չ����������xt(k,5)����', (xt(k,5),k=1,nxtemp)
	write(*,*)'��չ����������zt(k,1)����', (zt(k,1),k=1,nxtemp)
	write(*,*)'��չ����������zt(1,k)����', (zt(1,k),k=1,nytemp)
	pause
	!---------��������ǿ��������----------------------------------------
	do i=1,nxtemp
		do j=1,nytemp
			thetat(i,j)=xt(i,j)/(de/2.0)
			rt(i,j)=de/2.0-zt(i,j)
			xtnew1(i,j)=rt(i,j)*dcos(thetat(i,j))
			ytnew1(i,j)=rt(i,j)*dsin(thetat(i,j))
			ztnew1(i,j)=yt(i,j)
		end do
	end do
	!--------------------------------------------------------------------
	open(2,file='��������ǿ��������-3d.plt')  !��ʵ����һ��������������
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',2,'j=',nytemp,'k=',nxtemp
	do i=1,nxtemp
		do j=1,nytemp
			do k=1,2
				write(2,*) xtnew1(i,j),ytnew1(i,j),ztnew1(i,j)
			end do
		end do
	end do
	close(2)
	write(*,*)'��������ǿ��������-ztnew1(1,k)',(ztnew1(1,k),k=1,nytemp)
	pause
	!------------����Բ��------------------------------------------------
	do j=1,nxtemp
		do k=1,nytemp
			thetat(j,k)=xt(j,k)/(de/2.0)
			rt(j,k)=dm/2.0
			xtnew2(j,k)=rt(j,k)*dcos(thetat(j,k))
			ytnew2(j,k)=rt(j,k)*dsin(thetat(j,k))
			ztnew2(j,k)=yt(j,k)
		end do
	end do
	!--------------------------------------------------------------------
	open(2,file='����Բ��-3d.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nytemp,'j=',nxtemp,'k=',2
	do i=1,2
		do j=1,nxtemp
			do k=1,nytemp
				write(2,*) xtnew2(j,k),ytnew2(j,k),ztnew2(j,k)
			end do
		end do
	end do
	close(2)
	write(*,*)'����Բ��-ztnew2(1,k)',(ztnew2(1,k),k=1,nytemp)
	pause

	!--------------------------------------------------------------------
	allocate(xlast(301,301,301,24))
	allocate(ylast(301,301,301,24))   !��Ҫ�����ͷ��ڴ棬С����û�ˣ�������
	allocate(zlast(301,301,301,24))

	allocate(xlast1(301,301,301,22))
	allocate(ylast1(301,301,301,22))
	allocate(zlast1(301,301,301,22))
	!--------------------------------------------------------------------
	!***************************��һ����*********************************
	!�ɹ��������������
	!������������������κ���������С��
	jj1=int(nxtemp/4)/2+1
	jj2=jj1+int(nxtemp/4)
	jj3=jj2+int(nxtemp/4)
	jj4=jj3+int(nxtemp/4)
	!--------------------------------------------------------------------
	!����1
	m=1
	nxt(m)=9
	nyt(m)=int(nxtemp/4)+1
	nzt(m)=nytemp
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.4
	q=2.5
	call strech(n,p,q,s)
	do k=1,nzt(m)
		do j=1,nyt(m)
			if(j<=jj1) then
				xss(2,j)=xtnew1(jj4+j-1,k)
				yss(2,j)=ytnew1(jj4+j-1,k)
			else
				xss(2,j)=xtnew1(j-jj1+1,k)
				yss(2,j)=ytnew1(j-jj1+1,k)
			end if
			xss(1,j)=0.5*da
			yss(1,j)=-0.5*da+da*float(j-1)/float(nyt(m)-1)   !y����Ľڵ���δ����
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do

	do j= nyt(m)-10,nyt(m) !��Ϊ��10���������ȡֵ����11���ڵ�
		temporary2(3,j)=xss(1,j)
		temporary2(4,j)=yss(1,j)
	end do
	write(*,*) '����1ֱ�Ǳ���󼸸����x����temporary2(3,j)',(temporary2(3,j),j=nyt(m)-10,nyt(m))
	write(*,*) '����1ֱ�Ǳ���󼸸����y����temporary2(4,j)',(temporary2(4,j),j=nyt(m)-10,nyt(m))
	pause
	do j=nyt(m)-10,nyt(m)
		!����Щ���ݶ���ÿһά�ĵ�һ��Ԫ���ߣ�
		temporary3(1,j-int(nyt(m)/3)*2)=temporary2(3,j)
		temporary3(2,j-int(nyt(m)/3)*2)=temporary2(4,j)
	end do
	write(*,*) 'temporary3(1,j)=',(temporary3(1,j),j=1,int(nyt(1)/3))
	write(*,*) 'temporary3(2,j)=',(temporary3(2,j),j=1,int(nyt(1)/3))
	pause
	temporary(1)=xss(2,1)
	temporary(2)=yss(2,1)
	write(*,*)'����1�������߽翿��x���һ�����xֵtemporary(1)�ǣ�',temporary(1)
	write(*,*)'����1�������߽翿��x���һ�����yֵtemporary(2)�ǣ�',temporary(2)
	pause
	!--------------�������1����------------------------------------
	open(2,file='����1-3d.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				zlast(i,j,k,m)=yt(1,k)
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!����2
	m=2
	nxt(m)=nxt(1)
	nyt(m)=nyt(1)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.4
	q=2.5
	call strech(n,p,q,s)
	do k=1,nzt(m)
		do j=1,nyt(m)
			xss(2,j)=xtnew1(jj1+j-1,k)
			yss(2,j)=ytnew1(jj1+j-1,k)
			xss(1,j)= 0.50*da-da*float(j-1)/float(nyt(m)-1)
			yss(1,j)= 0.50*da
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	do j=1,nyt(m)
		temporary2(1,j)=xss(1,j)
		temporary2(2,j)=yss(1,j)
	end do
	write(*,*) '����2ֱ�Ǳ�x����',(temporary2(1,j),j=1,nyt(m))
	write(*,*) '����2ֱ�Ǳ�y����',(temporary2(2,j),j=1,nyt(m))
	pause
	!--------------------------------------------------------------------
	open(2,file='����2-3d.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!����3
	m=3
	nxt(m)=int(nyt(2)/3)
	nyt(m)=int(nyt(1)/3)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,2*n-1
		!��������
		!�ϱ߽�
		xs(1,2,i)=temporary2(1,i)
		ys(1,2,i)=temporary2(2,i) 
		!�±߽�
		xs(1,1,i)=xs(1,2,i)    
		ys(1,1,i)=temporary2(4,nyt(1)-10)
	end do

	temporary(13)=xs(1,1,1)
	temporary(15)=ys(1,1,1)
	temporary(14)=xs(1,1,nxt(m))
	
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do j=1,2*n-1
		!����������
		!��߽�
		xs(2,1,j)=xs(1,2,nxt(m))
		ys(2,1,j)=temporary3(2,j)    
		!�ұ߽�
		xs(2,2,j)=temporary3(1,j)   
		ys(2,2,j)=ys(2,1,j)
	end do
	!--------------------------------------------------------------------
	temporary(10)=xs(2,1,1)
	temporary(11)=ys(2,1,1)
	temporary(12)=ys(2,1,nyt(m))
	!   ���޲�ֵ��
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			xlast(i,j,1,m)=(1-sb(2,j))*xs(1,1,i)+sb(2,j)*xs(1,2,i)+(1-sb(1,i)) &
				&       *xs(2,1,j)+sb(1,i)*xs(2,2,j)-(sb(1,i)*sb(2,j)  &
				&       *xs(1,2,nxt(m))+sb(1,i)*(1-sb(2,j))*xs(1,1,nxt(m))+sb(2,j) &
				&       *(1-sb(1,i))*xs(1,2,1)+(1-sb(1,i)) &
				&       *(1-sb(2,j))*xs(1,1,1))

			ylast(i,j,1,m)=(1-sb(2,j))*ys(1,1,i)+sb(2,j)*ys(1,2,i)+(1-sb(1,i)) &
				&       *ys(2,1,j)+sb(1,i)*ys(2,2,j)-(sb(1,i)*sb(2,j) &
				&       *ys(1,2,nxt(m))+sb(1,i)*(1-sb(2,j))*ys(1,1,nxt(m))+sb(2,j) &
				&       *(1-sb(1,i))*ys(1,2,1)+(1-sb(1,i)) &
				&       *(1-sb(2,j))*ys(1,1,1))
		end do
	end do
	!--------------------------------------------------------------------
	open(2,file='����3-����-3ά.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=yt(1,k)
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!����4
	m=4
	nxt(m)=nyt(2)-nxt(3)+1
	nyt(m)=nyt(3)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do j=1,2*n-1
		xss(1,j)=0.0-dsin(pi/4.0)*da*dcos(pi/4.0)*sb(2,j)
		yss(1,j)=0.0+dsin(pi/4.0)*da*dsin(pi/4.0)*sb(2,j)
		xss(2,j)=temporary(10)
		yss(2,j)=temporary(11)+(temporary(12)-temporary(11))*sb(2,j)
	end do
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(1,i))*xss(1,j)+sb(1,i)*xss(2,j)
			ylast(i,j,1,m)=(1.0-sb(1,i))*yss(1,j)+sb(1,i)*yss(2,j)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����5
	m=5
	nxt(m)=nxt(3)
	nyt(m)=nxt(4)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,2*n-1
		xss(1,i)=0.0+dsin(pi/4.0)*da*dcos(pi/4.0)*sb(1,i)
		yss(1,i)=0.0-dsin(pi/4.0)*da*dsin(pi/4.0)*sb(1,i)
		xss(2,i)=temporary(14)+(temporary(13)-temporary(14))*sb(1,i)
		yss(2,i)=temporary(15)
	end do
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------

	open(2,file='��һ����������-2d.plt')
	WRITE(2,*) 'variables= "x","y","m"'
	!--------------------------------------------------------------------
	m=1
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=2
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=3
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=4
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=5
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!���������
	open(8,file='��һ����������-3d.plt')
	do m=1,5
		write(8,*) 'variables= "x","y","z","ms"'
		if(m==1) then
			write(8,*) 'zone t="domian01"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==2) then
			write(8,*) 'zone t="domian02"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==3) then
			write(8,*) 'zone t="domian03"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==4) then
			write(8,*) 'zone t="domian04"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==5) then
			write(8,*) 'zone t="domian05"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		do i=1,nxt(m)
			do j=1,nyt(m)
				do k=1,nzt(m)
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),m
				end do
			end do
		end do
	end do
	close(8)
	!-----------����һ���������������PLOT-3D��ʽ---------------------------
	nblocks=5
	do m=1,5
		imax(m)=nxt(m)
		jmax(m)=nyt(m)
		kmax(m)=nzt(m)
	end do
	open(8,file='plot3d-liutiyu-PCHE-3d.dat')
	write(8,*) nblocks
	write(8,*)(imax(m),jmax(m),kmax(m),m=1,5)
	do m=1,5
		write(8,*) &
			&    (((xlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
			&    (((ylast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
			&    (((zlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
	end do
	close(8)
	!--------------------------------------------------------------------
	!��������
	!����6
	m=6
	nxt(m)=7
	nyt(m)=nyt(1)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	!--------------------------------------------------------------------
	do k=1,nzt(m)
		do j=1,nyt(m)
			if(j<=jj1) then
				xss(2,j)=xtnew2(jj4+j-1,k)
				yss(2,j)=ytnew2(jj4+j-1,k)
				xss(1,j)=xtnew1(jj4+j-1,k)
				yss(1,j)=ytnew1(jj4+j-1,k)
			else
				xss(2,j)=xtnew2(j-jj1+1,k)
				yss(2,j)=ytnew2(j-jj1+1,k)
				xss(1,j)=xtnew1(j-jj1+1,k)
				yss(1,j)=ytnew1(j-jj1+1,k)
			end if
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	do j=1,(nyt(m)+1)/2
		temporary1(1,j)=xss(2,j)
		temporary1(2,j)=yss(2,j)
	end do
	write(*,*)'����6��Բ�߽��ǰһ��xֵtemporary1(1,j)�ǣ�',(temporary1(1,j),j=1,(nyt(m)+1)/2)
	write(*,*)'����6��Բ�߽��ǰһ��yֵtemporary1(2,j)�ǣ�',(temporary1(2,j),j=1,(nyt(m)+1)/2)
	pause
	!--------------------------------------------------------------------
	do j=(nyt(m)+1)/2,nyt(6)
		temporary1(10,j)=xss(2,j)
		temporary1(11,j)=yss(2,j)
	end do
	write(*,*)'����6��Բ�߽�ĺ�һ��xֵtemporary1(10,j)�ǣ�',(temporary1(10,j),j=(nyt(m)+1)/2,nyt(6))
	write(*,*)'����6��Բ�߽�ĺ�һ��yֵtemporary1(11,j)�ǣ�',(temporary1(11,j),j=(nyt(m)+1)/2,nyt(6))
	pause
	!--------------------------------------------------------------------
	open(2,file='����6-3ά.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!����7
	m=7
	nxt(m)=nxt(6)
	nyt(m)=nyt(2)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	!--------------------------------------------------------------------
	do k=1,nzt(m)
		do j=1,nyt(m)
			xss(2,j)=xtnew2(jj1+j-1,k)
			yss(2,j)=ytnew2(jj1+j-1,k)
			xss(1,j)=xtnew1(jj1+j-1,k)
			yss(1,j)=ytnew1(jj1+j-1,k)
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	do j=1,nyt(m)
		temporary1(3,j)=xss(2,j)
		temporary1(4,j)=yss(2,j)
	end do
	!--------------------------------------------------------------------
	temporary(3)=xss(1,nyt(m))
	temporary(4)=yss(1,nyt(m))
	write(*,*)'����7��Բ�߽��xֵtemporary1(3,j)�ǣ�',(temporary1(3,j),j=1,nyt(m))
	write(*,*)'����7��Բ�߽��yֵtemporary1(4,j)�ǣ�',(temporary1(4,j),j=1,nyt(m))
	write(*,*)'����7�������߽����һ����xֵtemporary(3)�ǣ�',temporary(3)
	write(*,*)'����7�������߽����һ����yֵtemporary(4)�ǣ�',temporary(4)
	pause
	!--------------------------------------------------------------------
	do j=1,(nyt(m)+1)/2
		!ȡ��Բ�߽�����һ��ֵ���������
		temporary3(5,j)=temporary1(3,nyt(7)+1-j)
		temporary3(6,j)=temporary1(4,nyt(7)+1-j)
	end do
	write(*,*)'temporary3(5,j)�ǣ�',(temporary3(5,j),j=1,(nyt(7)+1)/2)
	write(*,*)'temporary3(6,j)�ǣ�',(temporary3(6,j),j=1,(nyt(7)+1)/2)
	pause
	!--------------------------------------------------------------------
	do j=1,(nyt(m)+1)/2
		temporary1(12,j)=xss(2,j)
		temporary1(13,j)=yss(2,j)
	end do
	write(*,*)'����7��Բ�߽��ǰһ��xֵtemporary1(12,j)�ǣ�',(temporary1(12,j),j=1,(nyt(m)+1)/2)
	write(*,*)'����7��Բ�߽��ǰһ��yֵtemporary1(13,j)�ǣ�',(temporary1(13,j),j=1,(nyt(m)+1)/2)
	pause
	!--------------------------------------------------------------------
	open(2,file='����7-3ά.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!���������
	open(8,file='��һ���ֹ�����-3d.plt')
	do m=6,7
		write(8,*) 'variables= "x","y","z","ms"'
		if(m==6) then
			write(8,*) 'zone t="domian06"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==7) then
			write(8,*) 'zone t="domian07"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		do i=1,nxt(m)
			do j=1,nyt(m)
				do k=1,nzt(m)
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),m
				end do
			end do
		end do
	end do
	close(8)

	open(8,file='��һ����-3ά.plt')
	do m=1,7
		write(8,*) 'variables= "x","y","z","ms"'
		if(m==1) then
			write(8,*) 'zone t="domian01"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==2) then
			write(8,*) 'zone t="domian02"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==3) then
			write(8,*) 'zone t="domian03"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==4) then
			write(8,*) 'zone t="domian04"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==5) then
			write(8,*) 'zone t="domian05"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==6) then
			write(8,*) 'zone t="domian05"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==7) then
			write(8,*) 'zone t="domian05"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		do i=1,nxt(m)
			do j=1,nyt(m)
				do k=1,nzt(m)
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),m
				end do
			end do
		end do
	end do
	close(8)
	!***************�ڶ�����*********************************************
	!����8
	m=8
	nxt(m)=int(nxtemp/7)+1
	nyt(m)=(nyt(6)+1)/2
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	byta=dacos((dm+zd2)/((temp1)**2+(dm+zd2)**2)**0.5)
	do j=1,2*n-1
		xs(2,1,j)=temporary1(1,j)
		ys(2,1,j)=temporary1(2,j)
		xs(2,2,j)= (dm/2.0+temp1/2.0)*dcos(pi/4.0)+(temp1**2/4.0+dm**2/4.0+zd2**2/4.0+dm*zd2*0.5)**0.5*dcos(pi/4.0+byta)*sb(2,j)
		ys(2,2,j)=-(dm/2.0+temp1/2.0)*dsin(pi/4.0)+(temp1**2/4.0+dm**2/4.0+zd2**2/4.0+dm*zd2*0.5)**0.5*dsin(pi/4.0+byta)*sb(2,j)
	end do
	!--------------------------------------------------------------------
	temporary(16)=xs(2,2,nyt(m))
	temporary(17)=ys(2,2,nyt(m))
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(1,i))*xs(2,1,j)+sb(1,i)*xs(2,2,j)
			ylast(i,j,1,m)=(1.0-sb(1,i))*ys(2,1,j)+sb(1,i)*ys(2,2,j)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����9
	m=9
	nxt(m)=nxt(8)
	nyt(m)=(nyt(7)+1)/2
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do j=1,2*n-1
		xs(2,1,j)=temporary3(5,j)
		ys(2,1,j)=temporary3(6,j)
		xs(2,2,j)=-(dm/2.0+temp1/2.0)*dcos(pi/4.0)+(temp1**2/4.0+dm**2/4.0+zd2**2/4.0+dm*zd2*0.5)**0.5*dcos(pi/4.0-byta)*sb(2,j)
		ys(2,2,j)= (dm/2.0+temp1/2.0)*dsin(pi/4.0)+(temp1**2/4.0+dm**2/4.0+zd2**2/4.0+dm*zd2*0.5)**0.5*dsin(pi/4.0-byta)*sb(2,j)
	end do
	!--------------------------------------------------------------------
	temporary(18)=xs(2,2,nyt(m))
	temporary(19)=ys(2,2,nyt(m))
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(1,i))*xs(2,1,j)+sb(1,i)*xs(2,2,j)
			ylast(i,j,1,m)=(1.0-sb(1,i))*ys(2,1,j)+sb(1,i)*ys(2,2,j)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����10
	m=10
	nxt(m)=nyt(8)
	nyt(m)=11
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do j=1,2*n-1
		xss(1,j)= (dm/2.0+temp1/2.0)*dcos(pi/4.0)+(temp1/2.0)*dcos(pi/4.0)*sb(2,j)
		yss(1,j)=-(dm/2.0+temp1/2.0)*dsin(pi/4.0)-(temp1/2.0)*dsin(pi/4.0)*sb(2,j)
		xss(2,j)= temporary(16)+temp1*dcos(pi/4.0)*sb(2,j)
		yss(2,j)= temporary(17)-temp1*dsin(pi/4.0)*sb(2,j)
	end do
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(1,i))*xss(1,j)+sb(1,i)*xss(2,j)
			ylast(i,j,1,m)=(1.0-sb(1,i))*yss(1,j)+sb(1,i)*yss(2,j)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����11
	m=11
	nxt(m)=nyt(9)
	nyt(m)=nyt(10)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do j=1,2*n-1
		xss(1,j)=-(dm/2.0+temp1/2.0)*dcos(pi/4.0)-(temp1/2.0)*dcos(pi/4.0)*sb(2,j)
		yss(1,j)= (dm/2.0+temp1/2.0)*dsin(pi/4.0)+(temp1/2.0)*dsin(pi/4.0)*sb(2,j)
		xss(2,j)=temporary(18)-temp1*dcos(pi/4.0)*sb(2,j)
		yss(2,j)= temporary(19)+temp1*dsin(pi/4.0)*sb(2,j)
	end do
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(1,i))*xss(1,j)+sb(1,i)*xss(2,j)
			ylast(i,j,1,m)=(1.0-sb(1,i))*yss(1,j)+sb(1,i)*yss(2,j)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!**************��������**********************************************
	!����12
	m=12
	nxt(m)=7
	nyt(m)=nyt(10)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do

	!--------------------------------------------------------------------
	do j=1,2*n-1
		xss(1,j)=temporary(16)+temp1*dcos(pi/4.0)*sb(2,j)
		yss(1,j)= temporary(17)-temp1*dsin(pi/4.0)*sb(2,j)
		xss(2,j)=dm/2.0*dcos(pi/4.0)+(dm/2.0+zd2)*dcos(pi/4.0)+temp1*dcos(pi/4.0)*sb(2,j)
		yss(2,j)=-dm/2.0*dsin(pi/4.0)+(dm/2.0+zd2)*dsin(pi/4.0)-temp1*dsin(pi/4.0)*sb(2,j)
	end do
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(1,i))*xss(1,j)+sb(1,i)*xss(2,j)
			ylast(i,j,1,m)=(1.0-sb(1,i))*yss(1,j)+sb(1,i)*yss(2,j)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����13
	m=13
	nxt(m)=nxt(12)
	nyt(m)=nyt(11)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do

	!--------------------------------------------------------------------
	do j=1,2*n-1
		xss(1,j)=temporary(18)-temp1*dcos(pi/4.0)*sb(2,j)
		yss(1,j)=temporary(19)+temp1*dsin(pi/4.0)*sb(2,j)
		xss(2,j)=-dm/2.0*dcos(pi/4.0)+(dm/2.0+zd2)*dcos(pi/4.0)-temp1*dcos(pi/4.0)*sb(2,j)
		yss(2,j)=dm/2.0*dsin(pi/4.0)+(dm/2.0+zd2)*dsin(pi/4.0)+temp1*dsin(pi/4.0)*sb(2,j)
	end do
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(1,i))*xss(1,j)+sb(1,i)*xss(2,j)
			ylast(i,j,1,m)=(1.0-sb(1,i))*yss(1,j)+sb(1,i)*yss(2,j)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����14
	m=14
	nxt(m)=nxt(8)
	nyt(m)=((nyt(6)+1)/2+(nyt(7)+1)/2)-1
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do   j=1,nyt(m)
		sb(2,j)=s(j)
	end do
	!--------------------------------------------------------------------
	do j=1,(nyt(14)+1)/2
		temporary1(14,j)=temporary1(10,((nyt(14)+1)/2-1)+j)
		temporary1(15,j)=temporary1(11,((nyt(14)+1)/2-1)+j)
	end do
	write(*,*)'ǰһ��temporary1(14,j)',(temporary1(14,j),j=1,(nyt(14)+1)/2)
	write(*,*)'ǰһ��temporary1(15,j)',(temporary1(15,j),j=1,(nyt(14)+1)/2)
	pause
	!--------------------------------------------------------------------
	do j=1,(nyt(14)+1)/2
		temporary1(14,((nyt(14)+1)/2-1)+j)=temporary1(12,j)
		temporary1(15,((nyt(14)+1)/2-1)+j)=temporary1(13,j)
	end do
	write(*,*)'��һ��temporary1(14,j)',(temporary1(14,j),j=(nyt(14)+1)/2,nyt(14))
	write(*,*)'��һ��temporary1(15,j)',(temporary1(15,j),j=(nyt(14)+1)/2,nyt(14))
	pause
	!--------------------------------------------------------------------
	write(*,*)'������temporary1(14,j)',(temporary1(14,j),j=1,nyt(14))
	write(*,*)'������temporary1(15,j)',(temporary1(15,j),j=1,nyt(14))
	pause
	!--------------------------------------------------------------------
	do j=1,n
		xss(1,j)=temporary1(14,j)
		yss(1,j)= temporary1(15,j)
		xss(2,j)=dm/2.0*dcos(pi/4.0)+(dm/2.0+zd2/2.0)*dcos(pi/4.0)-dm*dcos(pi/4.0)*sb(2,j)
		yss(2,j)=-dm/2.0*dsin(pi/4.0)+(dm/2.0+zd2/2.0)*dsin(pi/4.0)+dm*dsin(pi/4.0)*sb(2,j)
	end do
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)

	do i=1,nxt(m)
		sb(1,i)=s(i)
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(1,i))*xss(1,j)+sb(1,i)*xss(2,j)
			ylast(i,j,1,m)=(1.0-sb(1,i))*yss(1,j)+sb(1,i)*yss(2,j)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����15
	m=15
	nxt(m)=nxt(12)
	nyt(m)=nyt(14)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,2*n-1
		xss(1,i)= dm/2.0*dcos(pi/4.0)+(dm/2.0+zd2/2.0)*dcos(pi/4.0)+zd2/2.0*dcos(pi/4.0)*sb(1,i)
		yss(1,i)=-dm/2.0*dsin(pi/4.0)+(dm/2.0+zd2/2.0)*dsin(pi/4.0)+zd2/2.0*dsin(pi/4.0)*sb(1,i)
		xss(2,i)=-dm/2.0*dcos(pi/4.0)+(dm/2.0+zd2/2.0)*dcos(pi/4.0)+zd2/2.0*dcos(pi/4.0)*sb(1,i)
		yss(2,i)= dm/2.0*dsin(pi/4.0)+(dm/2.0+zd2/2.0)*dsin(pi/4.0)+zd2/2.0*dsin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����16
	m=16
	nxt(m)=17
	nyt(m)=nyt(11)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,2*n-1
		xss(1,i)=-(dm/2.0+temp1/2.0)*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(1,i)=(dm/2.0+temp1/2.0)*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
		xss(2,i)=-(dm/2.0+temp1)*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(2,i)= (dm/2.0+temp1)*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����17
	m=17
	nxt(m)=nxt(16)
	nyt(m)=nxt(9)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,2*n-1
		xss(1,i)= -dm/2.0*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(1,i)=dm/2.0*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
		xss(2,i)=-(dm/2.0+temp1/2.0)*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(2,i)=(dm/2.0+temp1/2.0)*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����18
	m=18
	nyt(m)=nxt(7)
	nxt(m)= nxt(17)
	nzt(m)= nzt(1)
	n=nzt(m)
	do k=1,n
		yss(1,k)= -zt(1,k)
		zss(1,k)=  yt(1,k)
		yss(2,k)=(dm-de)/2.0
		zss(2,k)=yt(1,k)
	end do

	n=nyt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n
		sb(1,j)=s(j)
	end do

	do j=1,nyt(m)
		do k=1,nzt(m)
			ymm(j,k,m)=(1.0-sb(1,j))*yss(2,k)+sb(1,j)*yss(1,k)
			zmm(j,k,m)=(1.0-sb(1,j))*zss(2,k)+sb(1,j)*zss(1,k)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast1(i,j,k,m)=-temp3/(nxt(m)-1)*(i-1)
				ylast1(i,j,k,m)=ymm(j,k,m)
				zlast1(i,j,k,m)=zmm(j,k,m)
			end do
		end do
	end do
	!��ת+ƽ��
	!��ά
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=(xlast1(i,j,k,m)*dcos(-pi/4.0)+ylast1(i,j,k,m)*dsin(-pi/4.0))-(de/2.0)*dcos(pi/4.0)
				ylast(i,j,k,m)=(ylast1(i,j,k,m)*dcos(-pi/4.0)-xlast1(i,j,k,m)*dsin(-pi/4.0))+(de/2.0)*dsin(pi/4.0)
				zlast(i,j,k,m)=zlast1(i,j,k,m)
			end do
		end do
	end do

	open(2,file='train-2v.plt')
	write(2,*) 'variables= "y","z"'
	write(2,*) 'zone j=',nzt(m),'k=',nyt(m)
	do j=1,nyt(m)
		do k=1,nzt(m)
			write(2,*) ymm(j,k,m),zmm(j,k,m)
		end do
	end do
	close(2)

	open(2,file='train.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!����19
	m=19
	nyt(m)=nxt(1)
	nxt(m)= nxt(18)
	nzt(m)= nzt(1)
	n=nzt(m)
	do k=1,n
		yss(1,k)= -(de/2.0-dsin(pi/4.0)*da)
		zss(1,k)=  yt(1,k)
		yss(2,k)=-zt(1,k)
		zss(2,k)=yt(1,k)
	end do

	n=nyt(m)
	p=1.4
	q=2.5
	call strech(n,p,q,s)
	do j=1,n
		sb(1,j)=s(j)
	end do

	do j=1,nyt(m)
		do k=1,nzt(m)
			ymm(j,k,m)=(1.0-sb(1,j))*yss(1,k)+sb(1,j)*yss(2,k)
			zmm(j,k,m)=(1.0-sb(1,j))*zss(1,k)+sb(1,j)*zss(2,k)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast1(i,j,k,m)=-temp3/(nxt(m)-1)*(i-1)
				ylast1(i,j,k,m)=ymm(j,k,m)
				zlast1(i,j,k,m)=zmm(j,k,m)
			end do
		end do
	end do
	!��ת+ƽ��
	!��ά
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=(xlast1(i,j,k,m)*dcos(-pi/4.0)+ylast1(i,j,k,m)*dsin(-pi/4.0))-(de/2.0)*dcos(pi/4.0)
				ylast(i,j,k,m)=(ylast1(i,j,k,m)*dcos(-pi/4.0)-xlast1(i,j,k,m)*dsin(-pi/4.0))+(de/2.0)*dsin(pi/4.0)
				zlast(i,j,k,m)=zlast1(i,j,k,m)
			end do
		end do
	end do


	!--------------------------------------------------------------------
	!����20
	m=20
	nxt(m)=nxt(19)
	nyt(m)=nyt(4)+nxt(5)-1
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=int(nxt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do i=1,n-1
		s(n+i)=1.0-s(n-i)
	end do
	do   i=1,2*n-1
		if(i<=n)  then
			s(i)=s(i)/2.0
			sb(1,i)=s(i)
		else
			s(i)=0.5+s(i)/2.0
			sb(1,i)=s(i)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,2*n-1
		xss(1,i)= dsin(pi/4.0)*da*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(1,i)=-dsin(pi/4.0)*da*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
		xss(2,i)=-dsin(pi/4.0)*da*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(2,i)= dsin(pi/4.0)*da*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=int(nyt(m)/2)+1
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n-1
		s(n+j)=1.0-s(n-j)
	end do
	do  j=1,2*n-1
		if(j<=n)  then
			s(j)=s(j)/2.0
			sb(2,j)=s(j)
		else
			s(j)=0.5+s(j)/2.0
			sb(2,j)=s(j)
		end if
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����21
	m=21
	nyt(m)=nxt(2)
	nxt(m)= nxt(20)
	nzt(m)= nzt(1)
	n=nzt(m)
	do k=1,n
		yss(2,k)= -zt(1,k)
		zss(2,k)=yt(1,k)
		yss(1,k)=-(de/2.0-dsin(pi/4.0)*da)
		zss(1,k)=yt(1,k)
	end do

	n=nyt(m)
	p=1.4
	q=2.5
	call strech(n,p,q,s)
	do j=1,n
		sb(1,j)=s(j)
	end do
	do j=1,nyt(m)
		do k=1,nzt(m)
			ymm(j,k,m)=(1.0-sb(1,j))*yss(1,k)+sb(1,j)*yss(2,k)
			zmm(j,k,m)=(1.0-sb(1,j))*zss(1,k)+sb(1,j)*zss(2,k)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast1(i,j,k,m)=temp3/(nxt(m)-1)*(i-1)
				ylast1(i,j,k,m)=ymm(j,k,m)
				zlast1(i,j,k,m)=zmm(j,k,m)
			end do
		end do
	end do
	!��ת+ƽ��
	!��ά
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=(xlast1(i,j,k,m)*dcos(3*pi/4.0)+ylast1(i,j,k,m)*dsin(3*pi/4.0))+(de/2.0)*dcos(pi/4.0)
				ylast(i,j,k,m)=(ylast1(i,j,k,m)*dcos(3*pi/4.0)-xlast1(i,j,k,m)*dsin(3*pi/4.0))-(de/2.0)*dsin(pi/4.0)
				zlast(i,j,k,m)=zlast1(i,j,k,m)
			end do
		end do
	end do

	open(2,file='train-2v1.plt')
	write(2,*) 'variables= "y","z"'
	write(2,*) 'zone j=',nzt(m),'k=',nyt(m)
	do j=1,nyt(m)
		do k=1,nzt(m)
			write(2,*) ymm(j,k,m),zmm(j,k,m)
		end do
	end do
	close(2)

	open(2,file='train1.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!����22
	m=22
	nyt(m)=nxt(6)
	nxt(m)= nxt(21)
	nzt(m)= nzt(1)
	n=nzt(m)
	do k=1,n
		yss(1,k)= -zt(1,k)
		zss(1,k)=yt(1,k)
		yss(2,k)=(dm-de)/2.0
		zss(2,k)=yt(1,k)
	end do

	n=nyt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do j=1,n
		sb(1,j)=s(j)
	end do

	do j=1,nyt(m)
		do k=1,nzt(m)
			ymm(j,k,m)=(1.0-sb(1,j))*yss(1,k)+sb(1,j)*yss(2,k)
			zmm(j,k,m)=(1.0-sb(1,j))*zss(1,k)+sb(1,j)*zss(2,k)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast1(i,j,k,m)=temp3/(nxt(m)-1)*(i-1)
				ylast1(i,j,k,m)=ymm(j,k,m)
				zlast1(i,j,k,m)=zmm(j,k,m)
			end do
		end do
	end do
	!��ת+ƽ��
	!��ά
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=(xlast1(i,j,k,m)*dcos(3*pi/4.0)+ylast1(i,j,k,m)*dsin(3*pi/4.0))+(de/2.0)*dcos(pi/4.0)
				ylast(i,j,k,m)=(ylast1(i,j,k,m)*dcos(3*pi/4.0)-xlast1(i,j,k,m)*dsin(3*pi/4.0))-(de/2.0)*dsin(pi/4.0)
				zlast(i,j,k,m)=zlast1(i,j,k,m)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����23
	m=23
	nxt(m)=nxt(22)
	nyt(m)=nxt(8)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do   i=1,n
		sb(1,i)=s(i)
	end do
	!--------------------------------------------------------------------
	do i=1,n
		xss(1,i)= (dm/2.0+temp1/2.0)*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(1,i)=-(dm/2.0+temp1/2.0)*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
		xss(2,i)= dm/2.0*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(2,i)=-dm/2.0*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do   j=1,nyt(m)
		sb(2,j)=s(j)
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����24
	m=24
	nxt(m)=nxt(23)
	nyt(m)=nyt(10)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do   i=1,n
		sb(1,i)=s(i)
	end do
	!--------------------------------------------------------------------
	do i=1,n
		xss(1,i)= (dm/2.0+temp1)*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(1,i)=-(dm/2.0+temp1)*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
		xss(2,i)= (dm/2.0+temp1/2.0)*dcos(pi/4.0)-temp3*dcos(pi/4.0)*sb(1,i)
		yss(2,i)=-(dm/2.0+temp1/2.0)*dsin(pi/4.0)-temp3*dsin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	do   j=1,nyt(m)
		sb(2,j)=s(j)
	end do
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	do  i=1,nxt(m)
		do  j=1,nyt(m)
			do  k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!********************����PCHE��άͼ��*******************************
	open(2,file='PCHE-2d.plt')
	WRITE(2,*) 'variables= "x","y","m"'
	!--------------------------------------------------------------------
	m=1
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=2
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=3
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	m=4
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=5
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=6
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=7
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=8
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=9
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=10
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=11
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=12
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=13
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=14
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=15
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=16
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=17
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=18
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=19
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=20
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=21
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=22
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=23
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=24
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	close(2)
	!*********************�����άͼ��***********************************
	open(8,file='PCHE-3d.plt')
	do m=1,24
		write(8,*) 'variables= "x","y","z","ms"'
		if(m==1) then
			write(8,*) 'zone t="domian01"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==2) then
			write(8,*) 'zone t="domian02"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==3) then
			write(8,*) 'zone t="domian03"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==4) then
			write(8,*) 'zone t="domian04"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==5) then
			write(8,*) 'zone t="domian05"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==6) then
			write(8,*) 'zone t="domian06"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==7) then
			write(8,*) 'zone t="domian07"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==8) then
			write(8,*) 'zone t="domian08"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==9) then
			write(8,*) 'zone t="domian09"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==10) then
			write(8,*) 'zone t="domian10"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==11) then
			write(8,*) 'zone t="domian11"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==12) then
			write(8,*) 'zone t="domian12"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==13) then
			write(8,*) 'zone t="domian13"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==14) then
			write(8,*) 'zone t="domian14"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==15) then
			write(8,*) 'zone t="domian15"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==16) then
			write(8,*) 'zone t="domian16"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==17) then
			write(8,*) 'zone t="domian17"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==18) then
			write(8,*) 'zone t="domian18"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==19) then
			write(8,*) 'zone t="domian19"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==20) then
			write(8,*) 'zone t="domian20"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==21) then
			write(8,*) 'zone t="domian21"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==22) then
			write(8,*) 'zone t="domian22"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==23) then
			write(8,*) 'zone t="domian23"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==24) then
			write(8,*) 'zone t="domian24"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		do i=1,nxt(m)
			do j=1,nyt(m)
				do k=1,nzt(m)
					select case(m)
					case(1)
						mf=0
					case(2)
						mf=0
					case(3)
						mf=0
					case(4)
						mf=0
					case(5)
						mf=0
					case(6)
						mf=1
					case(7)
						mf=1
					case(8)
						mf=1
					case(9)
						mf=1
					case(10)
						mf=1
					case(11)
						mf=1
					case(12)
						mf=1
					case(13)
						mf=1
					case(14)
						mf=1
					case(15)
						mf=1
					case(16)
						mf=1
					case(17)
						mf=1
					case(18)
						mf=1
					case(19)
						mf=1
					case(20)
						mf=1
					case(21)
						mf=1
					case(22)
						mf=1
					case(23)
						mf=1
					case(24)
						mf=1
					end select
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),mf
				end do
			end do
		end do
	end do
	close(8)
	!*********************�����������άͼ��-Plot��ʽ**************************

	nblocks=19
	do m=6,24
		imax(m)=nxt(m)
		jmax(m)=nyt(m)
		kmax(m)=nzt(m)
	end do
	open(8,file='plot3d-gutiyu-PCHE-3d.dat')
	write(8,*) nblocks
	write(8,*)(imax(m),jmax(m),kmax(m),m=6,24)
	do m=6,24
		write(8,*) &
			&    (((xlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
			&    (((ylast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
			&    (((zlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
	end do
	close(8)

	end program main
	!**************���طֲ��������ӳ���**********************************
	subroutine strech(n,u,v,s)
	implicit none
	real(kind=8)  ::  s(301)
	integer  :: an,al,n,l
	real(kind=8)  :: deta,tqi,eta,dum,u,v
	an=n-1
	deta=1.0/float(an)
	tqi=1.0/tanh(v)
	do 1 l=1,n
		al=l-1
		eta=float(al)*deta
		dum=v*(1.0-eta)
		dum=1.0-tanh(dum)*tqi
		s(l)=u*eta+(1.0-u)*dum
1	continue
	return
	end subroutine strech
	!***********************OVER*********************************************