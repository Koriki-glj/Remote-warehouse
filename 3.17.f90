	!********************************START*******************************
	program main
	implicit none
	!********************************声明*******************************
	!边界和内部坐标
	real(kind=4):: xs(2,2,200),ys(2,2,200)
	real(kind=4),allocatable :: zlast(:,:,:,:),xlast(:,:,:,:),ylast(:,:,:,:)
	!分布函数
	real(kind=4) :: s(301),sb(2,301)
	!尺寸
	real(kind=4):: theta,theta3,theta4
	real(kind=4),parameter:: pi=4.0*atan(1.0)
	real(kind=4) :: zd2,temp1,temp3
	!网格正交性控制和其他循环变量
	real(kind=4) :: p,q
	integer :: n,j,m,k,i
	!Z方向的距离
	real(kind=4) :: lz
	real(kind=4) :: zjl(501,501)
	!-----------------金字塔————--------------------------------------
	!内部坐标
	real(kind=4):: yt(501,301),zt(501,301),xt(501,301)
	! 边界坐标
	real(kind=4):: xss(2,301),yss(2,301)
	!节点数
	integer :: nzt(31),nxt(31),nyt(31),imax(31),jmax(31),kmax(31)
	!PLOT-3D格式节点总数和颜色控制
	integer::  nblocks,mf
	!变量
	real(kind=4):: tempxp,tempyp
	!节点数
	integer:: nytemp,nxtemp,nztemp,nxra,nyrb,jj1,jj2,jj3,jj4
	integer:: nra,nrb,na,mb,nii,mjj
	!尺寸
	real(kind=4):: ra,rb,rh,alfa
	!临时一维数组
	real(kind=4):: temporary(20)
	!临时二维数组
	real(kind=4):: temporary1(8,201)
	!------------强化管外形--------------------------------------------------------
	real(kind=4):: ztnew1(501,301),xtnew1(501,301),ytnew1(501,301)
	real(kind=4):: rt(501,301),thetat(501,301)
	real(kind=4):: de,xL,da,yL
	!------------外置圆形---------------------------
	real(kind=4):: ztnew2(501,301),xtnew2(501,301),ytnew2(501,301)
	real(kind=4):: dm
	!********************************尺寸*******************************
	xL=4.2416
	yL=1.0531
	de=xL/pi
	da=0.6312*de
	nra=8
	nrb=3
	ra=xl/float(nra)
	rb=yl/float(nrb)
	alfa=atan(rb/ra)
	nztemp=31
	!以下参数独立定义
	rh=0.0338   !金字塔板的高度
	dm=1.50 !半圆的直径
	zd2=0.5
	temp1=0.5
	temp3=0.5
	lz=30
	!********************************金字塔建立***************************
	!金字塔粗糙元模型建立-双边界法
	!对X方向进行加密处理
	nxtemp=5
	n=int(nxtemp/2)+1  
	p=0.5285
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
	!写上下边界网格坐标
	do i=1,2*n-1
		xss(1,i)=(ra/2.0)*sb(1,i)
		yss(1,i)=0.0
		xss(2,i)=xss(1,i)
		yss(2,i)=rb/2.0
	end do
	!--------------------------------------------------------------------
	!对Y方向进行加密处理
	nytemp=11
	n=int(nytemp/2)+1
	p=0.5285
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
	!使用双边界法写内部网格
	do i=1,nxtemp
		do j=1,nytemp
			xt(i,j)=xss(1,i)*(1.0-sb(2,j))+sb(2,j)*xss(2,i)
			yt(i,j)=yss(1,i)*(1.0-sb(2,j))+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	!输出粗糙元网格
	open(2,file='金字塔-粗糙元网格-初始.plt')
	WRITE(2,*) 'variables= "x","y"'
	WRITE(2,*) 'zone i=',nytemp,'j=',nxtemp
	do  i=1,nxtemp
		do   j=1,nytemp
			write(2,*) xt(i,j),yt(i,j)
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!把粗糙元向x方向延长（网格疏密关系保持不变）,延长距离为ra/2.0
	do  i=1,nxtemp
		do   j=1,nytemp
			xt(nxtemp+i-1,j)=-xt(nxtemp-i+1,j)+ra
			yt(nxtemp+i-1,j)=yt(nxtemp-i+1,j)
		end do
	end do
	!--------------------------------------------------------------------
	!输出向x方向延长后的粗糙元网格
	open(2,file='金字塔-粗糙元网格+X.plt')
	write(2,*) 'variables= "x","y"'
	write(2,*) 'zone i=',nytemp,'j=',nxtemp
	do  i=1,nxtemp
		do   j=1,nytemp
			write(2,*) xt(i,j),yt(i,j)  !原先的部分
		end do
	end do

	write(2,*) 'zone i=',nytemp,'j=',nxtemp
	do  i=1,nxtemp
		do   j=1,nytemp
			write(2,*) xt(nxtemp+i-1,j),yt(nxtemp+i-1,j)  !延长的部分
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!把粗糙元向Y方向延长（网格疏密关系保持不变）,延长距离为rb/2.0
	nxra=2*nxtemp-1
	do  i=1,nxra
		do   j=1,nytemp
			xt(i,nytemp+j-1)=xt(i,nytemp-j+1)
			yt(i,nytemp+j-1)=-yt(i,nytemp-j+1)+rb
		end do
	end do
	!--------------------------------------------------------------------
	!输出向y方向延长后的粗糙元网格
	open(2,file='金字塔-粗糙元网格+Y.plt')
	WRITE(2,*) 'variables= "x","y"'
	WRITE(2,*) 'zone i=',nytemp,'j=',nxra
	do  i=1,nxra
		do   j=1,nytemp
			write(2,*) xt(i,j),yt(i,j)    !原先的部分
		end do
	end do
	!--------------------------------------------------------------------
	WRITE(2,*) 'zone i=',nytemp,'j=',nxra
	do  i=1,nxra
		do   j=1,nytemp
			write(2,*) xt(i,nytemp+j-1),yt(i,nytemp+j-1) !延长的部分
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!延长之后的整体粗糙元模型
	nyrb=2*nytemp-1
	open(2,file='金字塔-粗糙元+整体.plt')
	WRITE(2,*) 'variables= "x","y"'
	WRITE(2,*) 'zone i=',nyrb,'j=',nxra
	do  i=1,nxra
		do   j=1,nyrb
			write(2,*) xt(i,j),yt(i,j)
		end do
	end do
	close(2)
	!write(*,*)'yt',yt,这个为什么输出全是0；
	!pause
	write(*,*)'粗糙元-整体-yt(1,k)',(yt(1,k),k=1,nyrb)
	pause
	!--------------------------------------------------------------------
	!下面建立金字塔板模型
	do i=1,nxra
		do j=1,nyrb
			tempxp=abs(xt(i,j))  !返回参数的绝对值
			tempyp=abs(yt(i,j))
			!以下是粗糙元到金字塔板核心程序，但是看不懂
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
	open(2,file='金字塔板-初始.plt')
	WRITE(2,*) 'variables="x", "y","z"'
	WRITE(2,*) 'zone i=',2,'j=',nyrb,'k=',nxra   !这里的2会让y,z方向的数据每次重复2
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,2
				write(2,*) xt(i,j),yt(i,j),zt(i,j)
			end do
		end do
	end do
	close(2)
	write(*,*)'金字塔板-初始-yt(1,k)',(yt(1,k),k=1,nyrb)
	pause
	!--------------------------------------------------------------------
	!将金字塔板继续扩展
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
	!上述程序将整体金字塔板模型复制到一个新位置，其起点坐标为(nra-1)*ra，(nrb-1)*rb；
	open(2,file='金字塔板向右上角扩展.plt')
	WRITE(2,*) 'variables="x", "y","z"'
	WRITE(2,*) 'zone i=',1,'j=',nyrb,'k=',nxra
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,1
				write(2,*)  xt(nii+i,mjj+j),yt(nii+i,mjj+j),zt(nii+i,mjj+j)
			end do            !扩展后的部分
		end do
	end do
	!--------------------------------------------------------------------
	WRITE(2,*) 'zone i=',1,'j=',nyrb,'k=',nxra
	do i=1,nxra
		do  j=1,nyrb
			do   k=1,1
				write(2,*)   xt(i,j),yt(i,j),zt(i,j)
			end do            !扩展前的部分
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	!对X和Y方向的节点数重新赋值
	nxtemp=(nra-1)*(nxra-1)+nxra
	nytemp=(nrb-1)*(nyrb-1)+nyrb
	write(6,*)'扩展后的新节点数nxtemp和nytemp分别是：' ,nxtemp,nytemp
	pause
	!--------------------------------------------------------------------
	do i=1,nxtemp
		do j=1,nytemp
			zjl(i,j)=lz/(nytemp-1)*(j-1)
		end do
	end do
	
	write(*,*)'zjl(5,j)=',(zjl(5,j),j=1,nytemp)
	pause
	open(2,file='扩展后的金字塔板.plt')
	WRITE(2,*) 'variables="x", "y","z"'
	WRITE(2,*) 'zone i=',1,'j=',nytemp,'k=',nxtemp
	do i=1,nxtemp			 !扩展后的金字塔板长=nra*ra=xL,宽=nrb*rb=yL
		do  j=1,nytemp
			do   k=1,1
				write(2,*)xt(i,j),yt(i,j),zt(i,j)!zjl(j)
			end do
		end do
	end do
	close(2)
	write(*,*)'金字塔板-扩展-yt(1,k)',(yt(1,k),k=1,nytemp)
	pause
	!write(*,*)'扩展后金字塔板的yt(5,k)数据', (yt(5,k),k=1,nytemp)
	!write(*,*)'扩展后金字塔板的xt(k,5)数据', (xt(k,5),k=1,nxtemp)
	!write(*,*)'扩展后金字塔板的zt(k,5)数据', (zt(k,5),k=1,nxtemp)
	!write(*,*)'扩展后金字塔板的zt(5,k)数据', (zt(5,k),k=1,nytemp)
	!pause
	!---------金字塔型强化管外形----------------------------------------
	do i=1,nxtemp
		do j=1,nytemp
			thetat(i,j)=xt(i,j)/(de/2.0)
			rt(i,j)=de/2.0-zt(i,j)
			xtnew1(i,j)=rt(i,j)*cos(thetat(i,j))
			ytnew1(i,j)=rt(i,j)*sin(thetat(i,j))
			ztnew1(i,j)=zjl(i,j)
		end do
	end do
	!--------------------------------------------------------------------
	open(2,file='金字塔型强化管外形-3d.plt')  !其实做了一个流体域外轮廓
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
	write(*,*)'金字塔型强化管外形-ztnew1(1,k)',(ztnew1(1,k),k=1,nytemp)
	pause
	!------------外置圆形------------------------------------------------
	do j=1,nxtemp
		do k=1,nytemp
			thetat(j,k)=xt(j,k)/(de/2.0)
			rt(j,k)=dm/2.0
			xtnew2(j,k)=rt(j,k)*cos(thetat(j,k))
			ytnew2(j,k)=rt(j,k)*sin(thetat(j,k))
			ztnew2(j,k)=zjl(j,k)
		end do
	end do
	!--------------------------------------------------------------------
	open(2,file='外置圆形-3d.plt')
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
	write(*,*)'外置圆形-ztnew2(1,k)',(ztnew2(1,k),k=1,nytemp)
	pause
	!***************************第一部分*********************************
	!由固体域，流体域组成
	!流体域包括等腰三角形和两个流体小块
	!区块1
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
		xss(1,i)=0.0+da*sin(pi/4.0)*cos(pi/4.0)*sb(1,i)    !da是等腰三角形的腰长
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
	allocate(ylast(301,301,301,31))   !不要随意释放内存，小心数没了！！！！
	allocate(zlast(301,301,301,31))
	!--------------------------------------------------------------------
	do i=1,nxt(m)
		do j=1,nyt(m)
			xlast(i,j,1,m)=(1.0-sb(2,j))*xss(1,i)+sb(2,j)*xss(2,i)
			ylast(i,j,1,m)=(1.0-sb(2,j))*yss(1,i)+sb(2,j)*yss(2,i)
		end do
	end do
	!--------------------------------------------------------------------
	open(2,file='区块1-等腰三角形-2维.plt')
	write(2,*) 'variables= "x","y"'
	write(2,*) 'zone i=',nyt(m),'j=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			write(2,*) xlast(i,j,1,m),ylast(i,j,1,m)
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	open(2,file='区块1-等腰三角形-3维.plt')
	write(2,*) 'variables= "x","y","z"'
	write(2,*) 'zone i=',nzt(m),'j=',nyt(m),'k=',nxt(m)
	do i=1,nxt(m)
		do j=1,nyt(m)
			do k=1,nzt(m)
				xlast(i,j,k,m)=xlast(i,j,1,m)
				ylast(i,j,k,m)=ylast(i,j,1,m)
				zlast(i,j,k,m)=zjl(1,k)  !为何这里第三维拉伸用yt(1,k)的数据？,因为yt(1,k)一直到yt(481,k)的数据都是一样的
				write(2,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m)
			end do
		end do
	end do
	close(2)
	!--------------------------------------------------------------------
	open(2,file='yt(1,k)的数据.txt')
	write(2,*) (yt(1,k),k=1,nzt(m))  !隐含式循环
	close(2)
	!--------------------------------------------------------------------
	jj1=int(nxtemp/4)/2+1    !=61
	jj2=jj1+int(nxtemp/4)
	jj3=jj2+int(nxtemp/4)
	jj4=jj3+int(nxtemp/4)
	!--------------------------------------------------------------------
	!区块2
	m=2
	nxt(m)=11
	nyt(m)=nyt(1)
	nzt(m)=nzt(1)
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
			yss(1,j)=-0.5*da+da*float(j-1)/float(nyt(m)-1)   !y方向的节点数未加密
			do i=1,nxt(m)
				xlast(i,j,k,m)=(1.0-s(i))*xss(1,j)+s(i)*xss(2,j)
				ylast(i,j,k,m)=(1.0-s(i))*yss(1,j)+s(i)*yss(2,j)
			end do
		end do
	end do
	temporary(1)=xss(2,1)
	temporary(2)=yss(2,1)
	write(*,*)'区块2金字塔边界靠近x轴第一个点的x值temporary(1)是：',temporary(1)
	write(*,*)'区块2金字塔边界靠近x轴第一个点的y值temporary(2)是：',temporary(2)
	pause
	!--------------输出区块2网格------------------------------------
	open(2,file='区块2-3维.plt')
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
	open(2,file='zlast(1,1,k,1)的数据.txt')
	write(2,*) (zlast(1,1,k,1),k=1,nzt(m))
	close(2)
	!--------------------------------------------------------------------
	!区块3
	m=3
	nxt(m)=nxt(2)
	nyt(m)=nxt(1)
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
	open(2,file='区块3-3维.plt')
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
	!     OPEN(2,FILE='第一层高度.DAT')
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
	!输出流体域
	open(8,file='第一部分流体域-3维.plt')
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
	!-----------将第一部分流体域输出成PLOT-3D格式---------------------------
	!nblocks=3
	!do m=1,3
	!	imax(m)=nxt(m)
	!	jmax(m)=nyt(m)
	!	kmax(m)=nzt(m)
	!end do
	!open(8,file='plot3d-第一部分-fluid-3维.dat')
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
	!固体域建立
	!区块4
	m=4
	nxt(m)=7
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
	write(*,*)'区块4外圆边界的x值temporary1(1,j)是：',(temporary1(1,j),j=1,nyt(m))
	write(*,*)'区块4外圆边界的y值temporary1(2,j)是：',(temporary1(2,j),j=1,nyt(m))
	pause
	!--------------------------------------------------------------------
	open(2,file='区块4-3维.plt')
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
	!区块5
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
	write(*,*)'区块5外圆边界的x值temporary1(3,j)是：',(temporary1(3,j),j=1,nyt(m))
	write(*,*)'区块5外圆边界的y值temporary1(4,j)是：',(temporary1(4,j),j=1,nyt(m))
	write(*,*)'区块5金字塔边界靠近y轴的第一个点x值temporary(3)是：',temporary(3)
	write(*,*)'区块5金字塔边界靠近y轴的第一个点y值temporary(4)是：',temporary(4)
	pause
	!--------------------------------------------------------------------
	open(2,file='区块5-3维.plt')
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
	!输出固体域
	open(8,file='第一部分固体域-3维.plt')
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
	!将第一部分固体域输出成PLOT-3D格式
	!nblocks=2
	!do m=4,5
	!	imax(m)=nxt(m)
	!	jmax(m)=nyt(m)
	!	kmax(m)=nzt(m)
	!end do
	!open(8,file='plot3d-第一部分-solid-3维.dat')
	!write(8,*) nblocks
	!write(8,*)(imax(m),jmax(m),kmax(m),m=4,5)
	!do m=4,5
	!	write(8,*) &
	!		&    (((xlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
	!		&    (((ylast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
	!		&    (((zlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
	!end do
	!close(8)
	!输出第一部分三维图
	open(8,file='第一部分-3维.plt')
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
	!***************第二部分*********************************************
	!区块6
	m=6
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
		xs(2,2,j)= (dm/2.0)*cos(pi/4.0)+(dm/2.0+zd2)*cos(pi/4.0)+(dm/2.0)*cos(0.75*pi)*sb(2,j)
		ys(2,2,j)=-(dm/2.0)*sin(pi/4.0)+(dm/2.0+zd2)*sin(pi/4.0)+(dm/2.0)*sin(0.75*pi)*sb(2,j)
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
	!区块7
	m=7
	nxt(m)=nxt(6)
	nyt(m)=nyt(5)
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
		xs(2,2,j)=(dm/2.0+zd2)*cos(pi/4.0)+(dm/2.0)*cos(0.75*pi)*sb(2,j)
		ys(2,2,j)=(dm/2.0+zd2)*sin(pi/4.0)+(dm/2.0)*sin(0.75*pi)*sb(2,j)
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
	!区块8
	m=8
	nxt(m)=nxt(6)
	nyt(m)=17
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
		xss(1,j)= dm/2.0*cos(pi/4.0)+temp1*cos(pi/4.0)*sb(2,j)
		yss(1,j)=-dm/2.0*sin(pi/4.0)-temp1*sin(pi/4.0)*sb(2,j)
		xss(2,j)= dm/2.0*cos(pi/4.0)+(dm/2.0+zd2)*cos(pi/4.0)+temp1*cos(pi/4.0)*sb(2,j)
		yss(2,j)=-dm/2.0*sin(pi/4.0)+(dm/2.0+zd2)*sin(pi/4.0)-temp1*sin(pi/4.0)*sb(2,j)
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
	!区块9
	m=9
	nxt(m)=nxt(7)
	nyt(m)=nyt(8)
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
		xss(1,j)=-dm/2.0*cos(pi/4.0)-temp1*cos(pi/4.0)*sb(2,j)
		yss(1,j)= dm/2.0*sin(pi/4.0)+temp1*sin(pi/4.0)*sb(2,j)
		xss(2,j)=-dm/2.0*cos(pi/4.0)+(dm/2.0+zd2)*cos(pi/4.0)-temp1*cos(pi/4.0)*sb(2,j)
		yss(2,j)= dm/2.0*sin(pi/4.0)+(dm/2.0+zd2)*sin(pi/4.0)+temp1*sin(pi/4.0)*sb(2,j)
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
	!**************第五部分**********************************************
	!区块10
	m=10
	nxt(m)=13
	nyt(m)=nyt(9)
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
		xss(2,i)=-(dm/2.0+temp1)*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)= (dm/2.0+temp1)*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
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
	!区块11
	m=11
	nxt(m)=nxt(10)
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
	!区块12
	m=12
	nxt(m)=nxt(11)
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
	p=1.4
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
	!区块13
	m=13
	nxt(m)=nxt(12)
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

	!区块14
	m=14
	nxt(m)=nxt(13)
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
	p=1.4
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
	!区块15
	m=15
	nxt(m)=nxt(14)
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
	!区块16
	m=16
	nxt(m)=nxt(15)
	nyt(m)=nyt(8)
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
		xss(1,i)= (dm/2.0+temp1)*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(1,i)=-(dm/2.0+temp1)*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
		xss(2,i)= dm/2.0*cos(pi/4.0)-temp3*cos(pi/4.0)*sb(1,i)
		yss(2,i)=-dm/2.0*sin(pi/4.0)-temp3*sin(pi/4.0)*sb(1,i)
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
	!********************建立PCHE二维图像*******************************
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
	close(2)
	!*********************输出三维图像***********************************
	open(8,file='PCHE-3d.plt')
	do m=1,16
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
					end select
					write(8,*) xlast(i,j,k,m),ylast(i,j,k,m),zlast(i,j,k,m),mf
				end do
			end do
		end do
	end do
	close(8)
	!*********************输出三维图像-Plot格式**************************
	nblocks=16
	do m=1,16
		imax(m)=nxt(m)
		jmax(m)=nyt(m)
		kmax(m)=nzt(m)
	end do
	open(8,file='plot3d-PCHE-3d.dat')
	write(8,*) nblocks
	write(8,*)(imax(m),jmax(m),kmax(m),m=1,16)
	do m=1,16
		write(8,*) &
			&    (((xlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
			&    (((ylast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m)),&
			&    (((zlast(i,j,k,m),i=1,imax(m)),j=1,jmax(m)),k=1,kmax(m))
	end do
	close(8)
	end program main
	!**************返回分布函数的子程序**********************************
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