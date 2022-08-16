module array

    implicit none
    double precision, allocatable,dimension(:)::mi,th,ph,rwork,evl_s,n_i,s_x,s_y,s_z,n_i_up,n_i_dn
    integer,allocatable::x(:),y(:),z(:),tag(:),phase(:),b(:,:),b_temp(:),x_cl(:),y_cl(:)
    complex*16,allocatable:: h(:,:),bb(:,:),work(:),mat(:,:),h_(:,:),h1(:,:),aij(:,:),phi(:)
    double precision::nsum_t(1600),nn1(1600),ee1(1600),e_2(1600)

end module array

module global

    implicit none
    double precision::t,t2,ua,ub,uc,mu,m_sum,e_sum,u,esum_t,m1,ph1,th1,nsum,e,s_f,s_f1,s_f2,delta_t
    double precision::lambda,sx_i_sx_j,sy_i_sy_j,sz_i_sz_j,t_bi,t_ni_a,t_ni_b,t_ni,e_ni,e_o
    doubleprecision::n_ni,n_bi,n_o
    integer::nos,dim,unit_cells,nop,site,d_cl
    integer::i_1,i_2,i_3,k_x,k_y,b1,d,ie,sigma,a1,lx1,ly1,lx,ly,nx,ny,dkx,dky
    complex*16::ch_no,lambda1

endmodule global

program main
    use array
    use global
    implicit none
    read*,nop
    read*,d_cl
    read*,d
    read*,t
    read*,u
    site=d_cl**2
    dim=(nop)**site
    allocate(b((nop)**site,site),h(dim,dim))
    allocate(x_cl(site),y_cl(site),b_temp(site),evl_s(dim))
    allocate(x(d**2),y(d**2),h1(d**2,d**2),aij(d**2,d**2),phi(d**2))
    mu=0.4*u

    print*,'site in a cluster___________________________=',site
    print*,'nop on each site of cluster_________________=',nop
    print*,'dim of cluster H____________________________=',dim,'X',dim
    print*,'sites in square lattice for inter cluster___=',site
    print*,'dim of inter_cluster H______________________=',d**2,'X',d**2
    print*,'U___________________________________________=',u
    print*,'Mu__________________________________________=',mu

    
    call basis
    call cluster_matgen
    call diagonalization(h,1,dim)


    call inter_cluster_matgen
    call diagonalization(h1,2,d**2)
    call a_ij_cal
end program main

subroutine basis
    use array
    use global
    implicit none
    integer::i,j,k,l,i1,iy,ix,xi,xd,yi,yd,m,n,o,p,q
    
    i1=0
    do i=0,nop-1
        do j=0,nop-1
            do k=0,nop-1
                do l=0,nop-1
                    i1=i1+1                  
                    b(i1,4)=l
                    b(i1,3)=k
                    b(i1,2)=j
                    b(i1,1)=i
                    write(12,*) i1,b(i1,1),b(i1,2),b(i1,3),b(i1,4)                    
                enddo
            enddo
        enddo
    enddo


    ! i1=0
    ! do i=0,nop-1
    !     do j=0,nop-1
    !         do k=0,nop-1
    !             do l=0,nop-1
    !                 do m=0,nop-1
    !                     do n=0,nop-1
    !                         do o=0,nop-1
    !                             do p=0,nop-1
    !                                 do q=0,nop-1
    !                                     i1=i1+1                  
    !                                     b(i1,9)=q
    !                                     b(i1,8)=p
    !                                     b(i1,7)=o
    !                                     b(i1,6)=n
    !                                     b(i1,5)=m
    !                                     b(i1,4)=l
    !                                     b(i1,3)=k
    !                                     b(i1,2)=j
    !                                     b(i1,1)=i
    !                                     write(12,*) i1,b(i1,1),b(i1,2),b(i1,3),b(i1,4),b(i1,5),b(i1,6),b(i1,7)&
    !                                     &,b(i1,8),b(i1,9)
    !                                     flush(12)
    !                                 enddo
    !                             enddo
    !                         enddo
    !                     enddo
    !                 enddo
    !             enddo
    !         enddo
    !     enddo
    ! enddo
    !print*,b(1,2)
 
    i1=0
    do iy=1,d_cl
        do ix=1,d_cl
            i1=i1+1
            x_cl(i1)=ix
            y_cl(i1)=iy
        enddo
    enddo

    do l=1,d**2
        i=x_cl(l)
        j=y_cl(l)
        xi=1
        xd=-1
        yi=1                
        yd=-1 
  
        if (i.eq.1) xd=-i+(d_cl)
        if (i.eq.(d_cl)) xi=1-i
        if (j.eq.1) yd=-j+(d_cl)
        if (j.eq.(d_cl)) yi=1-j
        do k=1,d**2
            if(((x_cl(k).eq.(i+xi)).and.(y_cl(k).eq.j)))then
                write(500,*)l,k
            endif
            if(((x_cl(k).eq.i).and.(y_cl(k).eq.(j+yi))))then
                write(500,*)l,k
            endif              
            if(((x_cl(k).eq.(i+xd)).and.(y_cl(k).eq.j)))then
                write(500,*)l,k

            endif
            if(((x_cl(k).eq.i).and.(y_cl(k).eq.(j+yd))))then
                write(500,*)l,k

            endif 
        enddo
    enddo       

endsubroutine basis

subroutine cluster_matgen
    use array
    use global
    implicit none
    integer::i,j,i1,incre_n,decre_n,ix,xi,xd,yi,yd,l,k,j2,j1,b_1,b_2
    doubleprecision::n1,n2
    h=complex(0.0d0,0.0d0)

        do i_1=1,dim
            do i1=1,site
                b_temp(i1)=b(i_1,i1)
            enddo
            rewind(500)
            do i=1,site*4
                read(500,*)j1,j2
                !if(i_1==1) print*,j1,j2
                if(b_temp(j1).ne.0)then
                     !do i1=1,site
                     !enddo
                    b_1=b_temp(j1)-1
                    b_2=b_temp(j2)+1
                    n1=b_temp(j1)!*b_temp(j1)
                    n2=(b_temp(j2)+1)!*b_temp(j2)
                    b_temp(j1)=b_1
                    b_temp(j2)=b_2
                
                    do j=1,dim
                        i_2=0
                        do i1=1,site
                            if(b_temp(i1).eq.b(j,i1))then
                                i_2=i_2+1
                            endif
                        enddo
                        if(i_2.eq.site)then
                            h(i_1,j)=t*sqrt(n1*n2)
                            !_______________________________




                            !__________________________________
                            write(123,*) i_1,j
                            flush(123)
                        endif
                    enddo
                    b_temp(j1)=b_1+1
                    b_temp(j2)=b_2-1
                endif
            enddo
            b_temp=0
        enddo
  
    call matgen_u
    do i=1,dim
        do j=1,dim
            if(h(i,j).ne.conjg(h(j,i))) print*,i,j,h(i,j),h(j,i)
        enddo
    enddo
endsubroutine cluster_matgen

subroutine matgen_u
    use array
    use global
    implicit none
    integer::i,j,i1
    doubleprecision::n1,n2
    do i=1,dim
        n1=0
        n2=0
        do i1=1,site
            b_temp(i1)=b(i,i1)
            n1=n1+(b_temp(i1)-1)*b_temp(i1)
            n2=n2+b_temp(i1)
        enddo
        h(i,i)=u*n1/2-mu*n2
        b_temp=0
    enddo
endsubroutine matgen_u

subroutine diagonalization(h_temp,flag_diag,dim1)
    use array
    use global
    implicit none
    integer::lda,lwmax,info,lwork,flag_diag,i,dim1
    complex*16::h_temp(dim1,dim1)

    allocate(rwork(3*(dim1)-2),work(2*(dim1)-1))
    lda=(dim1)
    lwmax=(dim1)
    lwork=(2*(dim1)-1)
    evl_s=0.0d0

    if(flag_diag==1)then
       call zheev('n','u',dim1,h_temp,lda,evl_s,work,lwork,rwork,info)
    endif

    if(flag_diag==2)then
       call zheev('v','u',dim1,h_temp,lda,evl_s,work,lwork,rwork,info)
    endif

    if(info.ne.0)then
        print*,'algorithm failed'  
    endif 
   
    do i=1,dim1
        write(726,*) i, evl_s(i)
    enddo
    ! stop
    deallocate(rwork,work)  
endsubroutine diagonalization

subroutine dos
    use array
    use global
    implicit none
    double precision::pi,eta_,wmin,wmax,dw,w
    integer::nwp,j

    pi=4.0*atan(1.0)
    wmin=-30
    wmax=100
    nwp=2000
    dw=abs(wmax-wmin)/nwp
    w=wmin
    eta_=0.10d0
    do ie=1,nwp!no of interval bw bandwitdh of e	
        w=w+dw
        nsum=0.0d0

        do j=1,dim
           nsum=nsum+((eta_/pi)/((w-(evl_s(j)))**2+((eta_)**2)))
        enddo
   
        write(200,*) w,nsum/(dim)
        flush(200)
    enddo  
           
endsubroutine dos

subroutine inter_cluster_matgen
    use array
    use global
    implicit none
    integer::jx,jy,l,k,i,j,xi,xd,yd,yi
    t=-1
    j=0
    do jy=1,d
        do jx=1,d
            j=j+1
            x(j)=jx
            y(j)=jy
        enddo
    enddo

    do l=1,d**2
        if((mod(x(l),2).ne.0).and.((mod(y(l),2).ne.0)))then       
            i=x(l);j=y(l)          
            xi=1
            xd=-1
            yi=1                
            yd=-1 
            if (i.eq.1) xd=-i+d
            if (i.eq.d) xi=1-i
            if (j.eq.1) yd=-j+d
            if (j.eq.d) yi=1-j      
            do k=1,d**2           
                if(((x(k).eq.(i)).and.(y(k).eq.(j+yd))))then
                    write(501,*)l,k
                    h1(l,k)=t
                endif
                if(((x(k).eq.(i+xd)).and.(y(k).eq.(j))))then
                    write(501,*)l,k
                    h1(l,k)=t
                endif    
            enddo
        endif

        if((mod(x(l),2).ne.0).and.((mod(y(l),2).eq.0)))then       
            i=x(l);j=y(l)          
            xi=1
            xd=-1
            yi=1                
            yd=-1 
            if (i.eq.1) xd=-i+d
            if (i.eq.d) xi=1-i
            if (j.eq.1) yd=-j+d
            if (j.eq.d) yi=1-j      
            do k=1,d**2           
                if(((x(k).eq.(i)).and.(y(k).eq.(j+yi))))then
                    write(501,*)l,k
                    h1(l,k)=t
                endif
                if(((x(k).eq.(i+xd)).and.(y(k).eq.(j))))then
                    write(501,*)l,k
                    h1(l,k)=t
                endif    
            enddo
        endif

        if((mod(x(l),2).eq.0).and.((mod(y(l),2).ne.0)))then       
            i=x(l);j=y(l)          
            xi=1
            xd=-1
            yi=1                
            yd=-1 
            if (i.eq.1) xd=-i+d
            if (i.eq.d) xi=1-i
            if (j.eq.1) yd=-j+d
            if (j.eq.d) yi=1-j      
            do k=1,d**2           
                if(((x(k).eq.(i)).and.(y(k).eq.(j+yd))))then
                    write(501,*)l,k
                    h1(l,k)=t
                endif
                if(((x(k).eq.(i+xi)).and.(y(k).eq.(j))))then
                    write(501,*)l,k
                    h1(l,k)=t
                endif    
            enddo
        endif

        if((mod(x(l),2).eq.0).and.((mod(y(l),2).eq.0)))then       
            i=x(l);j=y(l)          
            xi=1
            xd=-1
            yi=1                
            yd=-1 
            if (i.eq.1) xd=-i+d
            if (i.eq.d) xi=1-i
            if (j.eq.1) yd=-j+d
            if (j.eq.d) yi=1-j      
            do k=1,d**2           
                if(((x(k).eq.(i)).and.(y(k).eq.(j+yi))))then
                    write(501,*)l,k
                    h1(l,k)=t
                endif
                if(((x(k).eq.(i+xi)).and.(y(k).eq.(j))))then
                    write(501,*)l,k
                    h1(l,k)=t
                endif    
            enddo
        endif

    enddo
endsubroutine inter_cluster_matgen

subroutine a_ij_cal
    use global
    use array
    implicit none
    integer::i,j,l
    complex*16::a_ij


    aij=complex(0.0d0,0.0d0)
    
    do i=1,d**2
        do j=1,d**2
            a_ij=complex(0.0d0,0.0d0)
            do l=1,d**2
                if(evl_s(l).lt.0)then
                    a_ij=a_ij+sqrt(abs(evl_s(l)))*conjg(h1(i,l))*h1(j,l)
                endif
            enddo
            aij(i,j)=a_ij
        enddo
    enddo

endsubroutine a_ij_cal

subroutine lambda_cal(p1)
    use global
    use array
    implicit none
    integer::i,j,l,p1
    complex*16::lambda_sum

    lambda_sum=0.0
    do j=1,d**2
        lambda_sum=lambda_sum+aij(p1,j)*phi(j)
    enddo
    lambda1=lambda_sum
endsubroutine lambda_cal

