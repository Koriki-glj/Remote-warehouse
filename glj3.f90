	!********************************START*******************************
	program main
	implicit none
	!********************************����*******************************
	!�߽���ڲ�����
	real(kind=4):: xs(2,2,200),ys(2,2,200)
	real(kind=4),allocatable :: zlast(:,:,:,:),xlast(:,:,:,:),ylast(:,:,:,:)
	!�ֲ�����
	real(kind=4) :: s(301),sb(2,301)
	!�ߴ�
	real(kind=4):: theta,theta3,theta4
	real(kind=4),parameter:: pi=4.0*atan(1.0)
	real(kind=4) :: zd2,temp1,temp2,temp3
	!���������Կ��ƺ�����ѭ������
	real(kind=4) :: p,q
	integer :: n,j,m,k,i
	!-----------------��������������--------------------------------------
	!�ڲ�����
	real(kind=4):: yt(501,201),zt(501,201),xt(501,201)
	! �߽�����
	real(kind=4):: xss(2,301),yss(2,301)
	!�ڵ���
	integer :: nzt(31),nxt(31),nyt(31),imax(31),jmax(31),kmax(31)
	!PLOT-3D��ʽ�ڵ���������ɫ����
	integer::  nblocks,mf
	!����
	real(kind=4):: tempxp,tempyp
	!�ڵ���
	integer:: nytemp,nxtemp,nxra,nyrb,jj1,jj2,jj3,jj4
	integer:: nra,nrb,na,mb,nii,mjj
	!�ߴ�
	real(kind=4):: ra,rb,rh,alfa
	!��ʱһά����
	real(kind=4):: temporary(20)
	!��ʱ��ά����
	real(kind=4):: temporary1(8,201)
	!------------ǿ��������--------------------------------------------------------
	real(kind=4):: ztnew1(501,201),xtnew1(501,201),ytnew1(501,201)
	real(kind=4):: rt(501,201),thetat(501,201)
	real(kind=4):: de,xL,da,yL
	!------------����Բ��---------------------------
	real(kind=4):: ztnew2(501,201),xtnew2(501,201),ytnew2(501,201)
	real(kind=4):: dm
	!********************************�ߴ�*******************************
	xL=207.24
	yL=51.81
	de=xL/pi
	da=0.6312*de
	nra=8
	nrb=3
	ra=xl/float(nra)
	rb=yl/float(nrb)
	rh=1.65
	alfa=atan(rb/ra)
	dm=70.0
	zd2=30.0
	temp1=20
	temp2=20
	temp3=20
	!********************************����������***************************
	!�������ֲ�Ԫģ�ͽ���-˫�߽編
	!��X������м��ܴ���
	nxtemp=31
	n=int(nxtemp/2)+1  !n=16
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
	nytemp=31
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
	WRITE(2,*) 'zone i=',2,'j=',nyrb,'k=',nxra
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,2
				write(2,*) xt(i,j),yt(i,j),zt(i,j)
			end do
		end do
	end do
	close(2)
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
	WRITE(2,*) 'zone i=',1,'j=',nyrb,'k=',nxra
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,1
				write(2,*)  xt(nii+i,mjj+j),yt(nii+i,mjj+j),zt(nii+i,mjj+j)
			end do            !��չ��Ĳ���
		end do
	end do
	!--------------------------------------------------------------------
	WRITE(2,*) 'zone i=',1,'j=',nyrb,'k=',nxra
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,1
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
	open(2,file='��չ��Ľ�������.plt')
	WRITE(2,*) 'variables="x", "y","z"'
	WRITE(2,*) 'zone i=',1,'j=',nytemp,'k=',nxtemp
	do i=1,nxtemp			 !��չ��Ľ������峤=nra*ra=xL,��=nrb*rb=yL
		do  j=1,nytemp
			do   k=1,1
				write(2,*)xt(i,j),yt(i,j),zt(i,j)
			end do
		end do
	end do
	close(2)
	!write(*,*)'��չ����������yt(5,k)����', (yt(5,k),k=1,nytemp)
	!write(*,*)'��չ����������zt(k,5)����', (zt(k,5),k=1,nxtemp)
	!write(*,*)'��չ����������zt(5,k)����', (zt(5,k),k=1,nytemp)
	!---------��������ǿ��������----------------------------------------
	do i=1,nxtemp
		do j=1,nytemp
			thetat(i,j)=xt(i,j)/(de/2.0)
			rt(i,j)=de/2.0-zt(i,j)
			xtnew1(i,j)=rt(i,j)*cos(thetat(i,j))
			ytnew1(i,j)=rt(i,j)*sin(thetat(i,j))
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
	!------------����Բ��------------------------------------------------
	do j=1,nxtemp
		do k=1,nytemp
			thetat(j,k)=xt(j,k)/(de/2.0)
			rt(j,k)=dm/2.0
			xtnew2(j,k)=rt(j,k)*cos(thetat(j,k))
			ytnew2(j,k)=rt(j,k)*sin(thetat(j,k))
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
	!***************************��һ����*********************************
	!�ɹ��������������
	!������������������κ���������С��
	!����1
	m=1
	nxt(1)=int(nxtemp/4)+1
	nyt(1)=nxt(1)
	nzt(1)=nytemp
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
		xss(1,i)=0.0+da*sin(pi/4.0)*cos(pi/4.0)*sb(1,i)
		yss(1,i)=0.0-da*sin(pi/4.0)*sin(pi/4.0)*sb(1,i)
		xss(2,i)=-0.5*da+da*sb(1,i)
		yss(2,i)=sin(pi/4.0)*da*sin(pi/4.0)
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
	allocate(xlast(301,301,301,31))
	allocate(ylast(301,301,301,31))   !��Ҫ�����ͷ��ڴ棬С����û�ˣ�������
	allocate(zlast(301,301,301,31))
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	open(2,file='����1-����������-2ά.plt')
	write(2,*) 'variables= "x","y"'
	write(2,*) 'zone i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m)
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	open(2,file='����1-����������-3ά.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=yt(1,k)  !Ϊ���������ά������yt(1,k)�����ݣ�,��Ϊyt(1,k)һֱ��yt(481,k)�����ݶ���һ����
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	open(2,file='yt(1,k)������.txt')
	write(2,*) (yt(1,k),k=1,nzt(m))
	close(2)
	!--------------------------------------------------------------------
	jj1=int(nxtemp/4)/2+1
	jj2=jj1+int(nxtemp/4)
	jj3=jj2+int(nxtemp/4)
	jj4=jj3+int(nxtemp/4)
	!--------------------------------------------------------------------
	!����2
	m=2
	nxt(m)=47
	nyt(m)=nyt(1)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n= nxt(m)
	p=1.6
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
			yss(1,j)=-0.5*da+da*float(j-1)/float(nyt(m)-1)
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	temporary(1)=xss(2,1)
	temporary(2)=yss(2,1)
	write(*,*)'����2�������߽翿��x���һ�����xֵtemporary(1)�ǣ�',temporary(1)
	write(*,*)'����2�������߽翿��x���һ�����yֵtemporary(2)�ǣ�',temporary(2)
	pause
	!--------------�������2����------------------------------------
	open(2,file='����2-3ά.plt')
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
	open(2,file='zlast(1,1,k,1)������.txt')
	write(2,*) (zlast(1,1,k,1),k=1,nzt(m))
	close(2)
	!--------------------------------------------------------------------
	!����3
	m=3
	nxt(m)=nxt(2)
	nyt(m)=nxt(1)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.6
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
	open(2,file='����3-3ά.plt')
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
	!     OPEN(2,FILE='��һ��߶�.DAT')
	!      WRITE(2,*) 'VARIABLES= "x","Y","Z"'
	!   WRITE(2,*) 'ZONE I=',NZT(M),'J=',NYT(M),'K=',NXT(M)
	!      do I=1,2
	!      DO J=1,NYT(2)
	!      DO K=1,NZT(2)
	!      WRITE(2,*) (Xlast(NXT(2),J,K,3)- Xlast(NXT(2)-1,J,K,3)),(ylast(NXT(2),J,K,3)- ylast(NXT(2)-1,J,K,3)),(zlast(NXT(2),J,K,3)- zlast(NXT(2)-1,J,K,3))
	!     END DO
	!   END DO
	!   end do
	!     CLOSE(2)
	!--------------------------------------------------------------------
	!���������
	open(8,file='��һ����������-3ά.plt')
	do m=1,3
		write(8,*) 'variables= "x","y","z"'
		if(m==1) then
			write(8,*) 'zone t="domian01"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==2) then
			write(8,*) 'zone t="domian02"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==3) then
			write(8,*) 'zone t="domian03"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		do i=1,nxt(m)
			do j=1,nyt(m)
				do k=1,nzt(m)
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
				end do
			end do
		end do
	end do
	close(8)
	!-----------����һ���������������PLOT-3D��ʽ---------------------------
	!nblocks=3
	!do m=1,3
	!	imax(m)=nxt(m)
	!	jmax(m)=nyt(m)
	!	kmax(m)=nzt(m)
	!end do
	!open(8,file='plot3d-��һ����-fluid-3ά.dat')
	!write(8,*) nblocks
	!write(8,*)(imax(m),jmax(m),kmax(m),m=1,3)
	!do m=1,3
	!	write(8,*) &
	!		&    (((xlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
	!		&    (((ylast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
	!		&    (((zlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
	!end do
	!close(8)
	!--------------------------------------------------------------------
	!��������
	!����4
	m=4
	nxt(m)=9
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
	do j=1,nyt(m)
		temporary1(1,j)=xss(2,j)
		temporary1(2,j)=yss(2,j)
	end do
	write(*,*)'����4��Բ�߽��xֵtemporary1(1,j)�ǣ�',(temporary1(1,j),j=1,nyt(m))
	write(*,*)'����4��Բ�߽��yֵtemporary1(2,j)�ǣ�',(temporary1(2,j),j=1,nyt(m))
	pause
	!--------------------------------------------------------------------
	open(2,file='����4-3ά.plt')
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
	!����5
	m=5
	nxt(m)=nxt(4)
	nyt(m)=nyt(3)
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
	write(*,*)'����5��Բ�߽��xֵtemporary1(3,j)�ǣ�',(temporary1(3,j),j=1,nyt(m))
	write(*,*)'����5��Բ�߽��yֵtemporary1(4,j)�ǣ�',(temporary1(4,j),j=1,nyt(m))
	write(*,*)'����5�������߽翿��y��ĵ�һ����xֵtemporary(3)�ǣ�',temporary(3)
	write(*,*)'����5�������߽翿��y��ĵ�һ����yֵtemporary(4)�ǣ�',temporary(4)
	pause
	!--------------------------------------------------------------------
	open(2,file='����5-3ά.plt')
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
	open(8,file='��һ���ֹ�����-3ά.plt')
	do m=4,5
		write(8,*) 'variables= "x","y","z"'
		if(m==4) then
			write(8,*) 'zone t="domian04"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==5) then
			write(8,*) 'zone t="domian05"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		do i=1,nxt(m)
			do j=1,nyt(m)
				do k=1,nzt(m)
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
				end do
			end do
		end do
	end do
	close(8)
	!����һ���ֹ����������PLOT-3D��ʽ
	!nblocks=2
	!do m=4,5
	!	imax(m)=nxt(m)
	!	jmax(m)=nyt(m)
	!	kmax(m)=nzt(m)
	!end do
	!open(8,file='plot3d-��һ����-solid-3ά.dat')
	!write(8,*) nblocks
	!write(8,*)(imax(m),jmax(m),kmax(m),m=4,5)
	!do m=4,5
	!	write(8,*) &
	!		&    (((xlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
	!		&    (((ylast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
	!		&    (((zlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
	!end do
	!close(8)
	!�����һ������άͼ
	open(8,file='��һ����-3ά.plt')
	do m=1,5
		write(8,*) 'variables= "x","y","z"'
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
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
				end do
			end do
		end do
	end do
	close(8)
	!***************�ڶ�����*********************************************
	!����6
	m=6
	nxt(m)=nxt(1)
	nyt(m)=nyt(1)
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
		xss(1,i)=0.0+da*sin(pi/4.0)*cos(pi/4.0)*sb(1,i)+(dm/2.0+zd2)*cos(pi/4.0)
		yss(1,i)=0.0-da*sin(pi/4.0)*sin(pi/4.0)*sb(1,i)+(dm/2.0+zd2)*sin(pi/4.0)
		xss(2,i)=-0.5*da+da*sb(1,i)+(dm/2.0+zd2)*cos(pi/4.0)
		yss(2,i)=sin(pi/4.0)*da*sin(pi/4.0)+(dm/2.0+zd2)*sin(pi/4.0)
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
	!����7
	m=7
	nxt(m)=nxt(2)
	nyt(m)=nyt(2)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n= nxt(m)
	p=1.6
	q=2.5
	call strech(n,p,q,s)
	do k=1,nzt(m)
		do j=1,nyt(m)
			if(j<=jj1) then
				xss(2,j)=xtnew1(jj4+j-1,k)+(dm/2.0+zd2)*cos(pi/4.0)
				yss(2,j)=ytnew1(jj4+j-1,k)+(dm/2.0+zd2)*sin(pi/4.0)
			else
				xss(2,j)=xtnew1(j-jj1+1,k)+(dm/2.0+zd2)*cos(pi/4.0)
				yss(2,j)=ytnew1(j-jj1+1,k)+(dm/2.0+zd2)*sin(pi/4.0)
			end if
			xss(1,j)=0.5*da+(dm/2.0+zd2)*cos(pi/4.0)
			yss(1,j)=-0.5*da+(dm/2.0+zd2)*sin(pi/4.0)+da*float(j-1)/float(nyt(m)-1)
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	temporary(5)=xss(1,1)
	temporary(6)=yss(1,1)
	temporary(7)=xss(2,1)
	temporary(8)=yss(2,1)
	write(*,*)'����7ֱ�߽��һ�����xֵtemporary(5)�ǣ�',temporary(5)
	write(*,*)'����7ֱ�߽��һ�����yֵtemporary(6)�ǣ�',temporary(6)
	write(*,*)'����7�������߽��һ�����xֵtemporary(7)�ǣ�',temporary(7)
	write(*,*)'����7�������߽��һ�����yֵtemporary(8)�ǣ�',temporary(8)
	pause
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����8
	m=8
	nxt(m)=nxt(3)
	nyt(m)=nyt(3)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.6
	q=2.5
	call strech(n,p,q,s)
	do k=1,nzt(m)
		do j=1,nyt(m)
			xss(2,j)=xtnew1(jj1+j-1,k)+(dm/2.0+zd2)*cos(pi/4.0)
			yss(2,j)=ytnew1(jj1+j-1,k)+(dm/2.0+zd2)*sin(pi/4.0)
			xss(1,j)=0.5*da+(dm/2.0+zd2)*cos(pi/4.0)-da*float(j-1)/float(nyt(m)-1)
			yss(1,j)=0.5*da+(dm/2.0+zd2)*sin(pi/4.0)
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	temporary(9)=xss(1,nyt(m))
	temporary(10)=yss(1,nyt(m))
	temporary(11)=xss(2,nyt(m))
	temporary(12)=yss(2,nyt(m))
	write(*,*)'����8ֱ�߽����һ�����xֵtemporary(9)�ǣ�',temporary(9)
	write(*,*)'����8ֱ�߽����һ�����yֵtemporary(10)�ǣ�',temporary(10)
	write(*,*)'����8�������߽����һ�����xֵtemporary(11)�ǣ�',temporary(11)
	write(*,*)'����8�������߽����һ�����yֵtemporary(12)�ǣ�',temporary(12)
	pause
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����9
	m=9
	nxt(m)=nxt(4)
	nyt(m)=nyt(4)
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
				xss(2,j)=xtnew2(jj4+j-1,k)+(dm/2.0+zd2)*cos(pi/4.0)
				yss(2,j)=ytnew2(jj4+j-1,k)+(dm/2.0+zd2)*sin(pi/4.0)
				xss(1,j)=xtnew1(jj4+j-1,k)+(dm/2.0+zd2)*cos(pi/4.0)
				yss(1,j)=ytnew1(jj4+j-1,k)+(dm/2.0+zd2)*sin(pi/4.0)
			else
				xss(2,j)=xtnew2(j-jj1+1,k)+(dm/2.0+zd2)*cos(pi/4.0)
				yss(2,j)=ytnew2(j-jj1+1,k)+(dm/2.0+zd2)*sin(pi/4.0)
				xss(1,j)=xtnew1(j-jj1+1,k)+(dm/2.0+zd2)*cos(pi/4.0)
				yss(1,j)=ytnew1(j-jj1+1,k)+(dm/2.0+zd2)*sin(pi/4.0)
			end if
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	temporary(13)=xss(2,1)
	temporary(14)=yss(2,1)
	do j=1,nyt(m)
		temporary1(5,j)=xss(2,j)
		temporary1(6,j)=yss(2,j)
	end do
	!--------------------------------------------------------------------
	write(*,*)'����9��Բ�߽��һ�����xֵtemporary(13)�ǣ�',temporary(13)
	write(*,*)'����9��Բ�߽��һ�����yֵtemporary(14)�ǣ�',temporary(14)
	write(*,*)'����9��Բ�߽��xֵtemporary1(5,j)�ǣ�',(temporary1(5,j),j=1,nyt(m))
	write(*,*)'����9��Բ�߽��yֵtemporary1(6,j)�ǣ�',(temporary1(6,j),j=1,nyt(m))
	pause
	!--------------------------------------------------------------------

	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	!����10
	m=10
	nxt(m)=nxt(5)
	nyt(m)=nyt(5)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	!--------------------------------------------------------------------
	do k=1,nzt(m)
		do j=1,nyt(m)
			xss(2,j)=xtnew2(jj1+j-1,k)+(dm/2.0+zd2)*cos(pi/4.0)
			yss(2,j)=ytnew2(jj1+j-1,k)+(dm/2.0+zd2)*sin(pi/4.0)
			xss(1,j)=xtnew1(jj1+j-1,k)+(dm/2.0+zd2)*cos(pi/4.0)
			yss(1,j)=ytnew1(jj1+j-1,k)+(dm/2.0+zd2)*sin(pi/4.0)
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	!--------------------------------------------------------------------
	temporary(15)=xss(2,nyt(m))
	temporary(16)=yss(2,nyt(m))
	do j=1,nyt(m)
		temporary1(7,j)=xss(2,j)
		temporary1(8,j)=yss(2,j)
	end do

	!--------------------------------------------------------------------
	write(*,*)'����10��Բ�߽翿��y���һ�����xֵtemporary(15)�ǣ�',temporary(15)
	write(*,*)'����10��Բ�߽翿��y���һ�����yֵtemporary(16)�ǣ�',temporary(16)
	write(*,*)'����10��Բ�߽��xֵtemporary1(7,j)�ǣ�',(temporary1(7,j),j=1,nyt(m))
	write(*,*)'����10��Բ�߽��yֵtemporary1(8,j)�ǣ�',(temporary1(8,j),j=1,nyt(m))
	pause
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				zlast(i,j,k,m)=zlast(1,1,k,1)
			end do
		end do
	end do
	!********************************��������*******************************
	!����11
	m=11
	nxt(m)=int(nxtemp/4)+1
	nyt(m)=nyt(9)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	!--------------------------------------------------------------------
	do k=1,nzt(m)
		do j=1,nyt(m)
			xs(2,1,j)=temporary1(5,j)
			ys(2,1,j)=temporary1(6,j)
			xs(2,2,j)= dm/2.0*cos(pi/4.0)+(dm+zd2+temp2)*cos(pi/4.0)+(dm/2.0)*cos(0.75*pi)*sb(2,j)
			ys(2,2,j)=-dm/2.0*sin(pi/4.0)+(dm+zd2+temp2)*sin(pi/4.0)+(dm/2.0)*sin(0.75*pi)*sb(2,j)
			do i=1,nxt(m)
				xlast(i,j,1,m)=(1.0-s(i))*xs(2,1,j)+s(i)*xs(2,2,j)
				ylast(i,j,1,m)=(1.0-s(i))*ys(2,1,j)+s(i)*ys(2,2,j)
			end do
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
	!����12
	m=12
	nxt(m)=nxt(11)
	nyt(m)=nyt(10)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nxt(m)
	p=1.0
	q=2.0
	call strech(n,p,q,s)
	!--------------------------------------------------------------------
	do k=1,nzt(m)
		do j=1,nyt(m)
			xs(2,1,j)=temporary1(7,j)
			ys(2,1,j)=temporary1(8,j)
			xs(2,2,j)=dm/2.0*cos(pi/4.0)+(dm/2.0+zd2+temp2)*cos(pi/4.0)+(dm/2.0)*cos(0.75*pi)*sb(2,j)
			ys(2,2,j)=dm/2.0*sin(pi/4.0)+(dm/2.0+zd2+temp2)*sin(pi/4.0)+(dm/2.0)*sin(0.75*pi)*sb(2,j)
			do i=1,nxt(m)
				xlast(i,j,1,m)=(1.0-s(i))*xs(2,1,j)+s(i)*xs(2,2,j)
				ylast(i,j,1,m)=(1.0-s(i))*ys(2,1,j)+s(i)*ys(2,2,j)
			end do
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
	nxt(m)=nxt(11)
	nyt(m)=int(nxtemp/8)+1
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
		xss(1,j)= dm/2.0*cos(pi/4.0)+(dm/2.0+zd2)*cos(pi/4.0)+temp1*cos(pi/4.0)*sb(2,j)
		yss(1,j)=-dm/2.0*sin(pi/4.0)+(dm/2.0+zd2)*sin(pi/4.0)-temp1*sin(pi/4.0)*sb(2,j)
		xss(2,j)= dm/2.0*cos(pi/4.0)+(dm+zd2+temp2)*cos(pi/4.0)+temp1*cos(pi/4.0)*sb(2,j)
		yss(2,j)=-dm/2.0*sin(pi/4.0)+(dm+zd2+temp2)*sin(pi/4.0)-temp1*sin(pi/4.0)*sb(2,j)
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
	nxt(m)=nxt(12)
	nyt(m)=nyt(13)
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
		xss(1,j)=-dm/2.0*cos(pi/4.0)+(dm/2.0+zd2)*cos(pi/4.0)-temp1*cos(pi/4.0)*sb(2,j)
		yss(1,j)= dm/2.0*sin(pi/4.0)+(dm/2.0+zd2)*sin(pi/4.0)+temp1*sin(pi/4.0)*sb(2,j)
		xss(2,j)=-dm/2.0*cos(pi/4.0)+(dm/2.0+zd2+dm/2.0+temp2)*cos(pi/4.0)-temp1*cos(pi/4.0)*sb(2,j)
		yss(2,j)= dm/2.0*sin(pi/4.0)+(dm/2.0+zd2+dm/2.0+temp2)*sin(pi/4.0)+temp1*sin(pi/4.0)*sb(2,j)
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
	!**************���Ĳ���**********************************************
	!����15
	m=15
	nxt(m)=int(nxtemp/4)+1
	nyt(m)=nyt(4)
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
		xs(2,1,j)=temporary1(1,j)
		ys(2,1,j)=temporary1(2,j)
		xs(2,2,j)=temporary(5)+sin(pi/4.0)*da*cos(0.75*pi)*sb(2,j)
		ys(2,2,j)=temporary(6)+sin(pi/4.0)*da*sin(0.75*pi)*sb(2,j)
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
	!����16
	m=16
	nxt(m)=nxt(15)
	nyt(m)=nyt(6)
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
		xs(2,1,j)=temporary1(3,j)
		ys(2,1,j)=temporary1(4,j)
		xs(2,2,j)=(dm/2.0+zd2)*cos(pi/4.0)+sin(pi/4.0)*da*cos(0.75*pi)*sb(2,j)
		ys(2,2,j)=(dm/2.0+zd2)*sin(pi/4.0)+sin(pi/4.0)*da*sin(0.75*pi)*sb(2,j)
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
	!����17
	m=17
	nxt(m)=nxt(15)
	nyt(m)=nxt(7)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.6
	q=2.5
	call strech(n,p,q,s)
	do   j=1,nyt(m)
		sb(2,j)=s(j)
	end do
	!--------------------------------------------------------------------
	do j=1,2*n-1
		xss(2,j)= dm/2.0*cos(pi/4.0)+temp1*cos(pi/4.0)*sb(2,j)
		yss(2,j)=-dm/2.0*sin(pi/4.0)-temp1*sin(pi/4.0)*sb(2,j)
		xss(1,j)=temporary(5)+(temporary(7)-temporary(5))*sb(2,j)
		yss(1,j)=temporary(6)+(temporary(8)-temporary(6))*sb(2,j)
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
	!����18
	m=18
	nxt(m)=nxt(17)
	nyt(m)=nxt(9)
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
	do j=1,2*n-1
		xss(1,j)= (dm/2.0+temp1)*cos(pi/4.0)+dm/8.0*cos(pi/4.0)*sb(2,j)
		yss(1,j)=-(dm/2.0+temp1)*sin(pi/4.0)+dm/8.0*sin(pi/4.0)*sb(2,j)
		xss(2,j)=temporary(7)+(temporary(13)-temporary(7))*sb(2,j)
		yss(2,j)=temporary(8)+(temporary(14)-temporary(8))*sb(2,j)
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
	temporary(17)=xss(1,nyt(m))
	temporary(18)=yss(1,nyt(m))
	write(*,*)'����18��߽����һ�����xֵtemporary(17)�ǣ�',temporary(17)
	write(*,*)'����18��߽����һ�����yֵtemporary(18)�ǣ�',temporary(18)
	pause
	!--------------------------------------------------------------------
	!����19
	m=19
	nxt(m)=nxt(18)
	nyt(m)=nyt(13)
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
		xss(1,j)=temporary(17)+(3.0*dm/16.0+zd2/2.0)*cos(pi/4.0)*sb(2,j)
		yss(1,j)=temporary(18)+(3.0*dm/16.0+zd2/2.0)*sin(pi/4.0)*sb(2,j)
		xss(2,j)=temporary(13)+temp1*cos(pi/4.0)*sb(2,j)
		yss(2,j)=temporary(14)-temp1*sin(pi/4.0)*sb(2,j)
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
	!����20
	m=20
	nxt(m)=nxt(16)
	nyt(m)=nxt(8)
	nzt(m)=nzt(1)
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.6
	q=2.5
	call strech(n,p,q,s)
	do   j=1,nyt(m)
		sb(2,j)=s(j)
	end do
	!--------------------------------------------------------------------
	do j=1,2*n-1
		xss(1,j)=-dm/2.0*cos(pi/4.0)-temp1*cos(pi/4.0)*sb(2,j)
		yss(1,j)= dm/2.0*sin(pi/4.0)+temp1*sin(pi/4.0)*sb(2,j)
		xss(2,j)=temporary(9)+(temporary(11)-temporary(9))*sb(2,j)
		yss(2,j)=temporary(10)+(temporary(12)-temporary(10))*sb(2,j)
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
	!����21
	m=21
	nxt(m)=nxt(20)
	nyt(m)=nxt(10)
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
	do j=1,2*n-1
		xss(1,j)=-(dm/2.0+temp1)*cos(pi/4.0)+dm/8.0*cos(pi/4.0)*sb(2,j)
		yss(1,j)= (dm/2.0+temp1)*sin(pi/4.0)+dm/8.0*sin(pi/4.0)*sb(2,j)
		xss(2,j)=temporary(11)+(temporary(15)-temporary(11))*sb(2,j)
		yss(2,j)=temporary(12)+(temporary(16)-temporary(12))*sb(2,j)
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
	temporary(19)=xss(1,nyt(m))
	temporary(20)=yss(1,nyt(m))
	write(*,*)'����21��߽����һ�����xֵtemporary(19)�ǣ�',temporary(19)
	write(*,*)'����21��߽����һ�����yֵtemporary(20)�ǣ�',temporary(20)
	pause
	!--------------------------------------------------------------------
	!����22
	m=22
	nxt(m)=nxt(21)
	nyt(m)=nyt(14)
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
		xss(1,j)=temporary(19)+(3.0*dm/16.0+zd2/2.0)*cos(pi/4.0)*sb(2,j)
		yss(1,j)=temporary(20)+(3.0*dm/16.0+zd2/2.0)*sin(pi/4.0)*sb(2,j)
		xss(2,j)=temporary(15)-temp1*cos(pi/4.0)*sb(2,j)
		yss(2,j)=temporary(16)+temp1*sin(pi/4.0)*sb(2,j)
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
	!**************���岿��**********************************************
	!����23
	m=23
	nxt(m)=45
	nyt(m)=nyt(20)
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
		xss(1,i)=-dm/2.0*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(1,i)= dm/2.0*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
		xss(2,i)=(-dm/2.0*cos(pi/4.0)-temp1*cos(pi/4.0))-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)=( dm/2.0*sin(pi/4.0)+temp1*cos(pi/4.0))-temp3*sin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.6
	q=2.5
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
	nyt(m)=nxt(5)
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
		xss(1,i)=temporary(3)-temp3*cos(pi/4.0)*sb(1,i)
		yss(1,i)=temporary(4)-temp3*sin(pi/4.0)*sb(1,i)
		xss(2,i)=-dm/2.0*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)= dm/2.0*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
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
	!����25
	m=25
	nxt(m)=nxt(24)
	nyt(m)=nxt(3)
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
		xss(1,i)=-sin(pi/4.0)*da*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(1,i)= sin(pi/4.0)*da*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
		xss(2,i)=temporary(3)-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)=temporary(4)-temp3*sin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.6
	q=2.5
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
	!����26
	m=26
	nxt(m)=nxt(25)
	nyt(m)=nxt(1)+nyt(1)-1
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
		xss(1,i)= sin(pi/4.0)*da*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(1,i)=-sin(pi/4.0)*da*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
		xss(2,i)=-sin(pi/4.0)*da*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)= sin(pi/4.0)*da*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
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

	!����27
	m=27
	nxt(m)=nxt(26)
	nyt(m)=nxt(2)
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
		xss(2,i)=temporary(1)-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)=temporary(2)-temp3*sin(pi/4.0)*sb(1,i)
		xss(1,i)= sin(pi/4.0)*da*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(1,i)=-sin(pi/4.0)*da*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.6
	q=2.5
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
	!����28
	m=28
	nxt(m)=nxt(27)
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
		xss(1,i)= dm/2.0*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(1,i)=-dm/2.0*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
		xss(2,i)=temporary(1)-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)=temporary(2)-temp3*sin(pi/4.0)*sb(1,i)
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
	!����29
	m=29
	nxt(m)=nxt(28)
	nyt(m)=nyt(17)
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
		xss(2,i)= (dm/2.0+temp1)*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)=-(dm/2.0+temp1)*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
		xss(1,i)= dm/2.0*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(1,i)=-dm/2.0*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
	end do
	!--------------------------------------------------------------------
	n=nyt(m)
	p=1.6
	q=2.5
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
	!--------------------------------------------------------------------
	m=4
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=5
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=6
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=7
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=8
	write(2,*) 'zone t="fluid"','i=',nyt(m),'j=',nxt(m)
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
	!--------------------------------------------------------------------
	m=25
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=26
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=27
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=28
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	!--------------------------------------------------------------------
	m=29
	write(2,*) 'zone t="solid"','i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m),m
		end do
	end do
	close(2)
	!*********************�����άͼ��***********************************
	open(8,file='PCHE-3d.plt')
	do m=1,29
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
		if(m==25) then
			write(8,*) 'zone t="domian25"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==26) then
			write(8,*) 'zone t="domian26"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==27) then
			write(8,*) 'zone t="domian27"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==28) then
			write(8,*) 'zone t="domian28"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
		end if
		if(m==29) then
			write(8,*) 'zone t="domian29"','i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
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
						mf=1
					case(5)
						mf=1
					case(6)
						mf=0
					case(7)
						mf=0
					case(8)
						mf=0
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
					case(25)
						mf=1
					case(26)
						mf=1
					case(27)
						mf=1
					case(28)
						mf=1
					case(29)
						mf=1
					end select
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),mf
				end do
			end do
		end do
	end do
	close(8)
	end program main
	!**************���طֲ��������ӳ���**********************************
	subroutine strech(n,u,v,s)
	implicit none
	real(kind=4)  ::  s(301)
	integer  :: an,al,n,l
	real(kind=4)  :: deta,tqi,eta,dum,u,v
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