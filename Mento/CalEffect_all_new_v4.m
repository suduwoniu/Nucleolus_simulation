%This function is for calculation of binding ability of each processing factor

function[randomFactor_1]=CalEffect_all_new_v4(randomFactor,inner,n1,n2,n3,range)
s=size(randomFactor);
%position: record the position of each processing factor
position=zeros(n1,n2,n3);
for n=1:s(1)
    position(randomFactor(n,2),randomFactor(n,3),randomFactor(n,4))=1;
end

%calculation of the number of processing factors in each 3*3*3 space, this
%is the processing density of each 3*3*3 space
dense=zeros(n1,n2,n3);
for xx=1:n1
    for yy=1:n2
        for zz=1:n3
            count0=0;
             for ii=-1:1
                  for jj=-1:1
                      for kk=-1:1
                           if xx+ii>0 && xx+ii<=n1 && yy+jj>0 && yy+jj<=n2 && zz+kk>0 && zz+kk<=n3
                               if position(xx+ii,yy+jj,zz+kk)==1
                                 count0=count0+1;
                               end
                           end
                      end
                  end
             end
             dense(xx,yy,zz)=count0;
        end
    end
end
dense_1=dense;
for xx=1:n1
    for yy=1:n2
        for zz=1:n3
            if dense(xx,yy,zz)<15
                cn=0;
                 for ii=-1:1
                  for jj=-1:1
                      for kk=-1:1
                           if xx+ii>0 && xx+ii<=n1 && yy+jj>0 && yy+jj<=n2 && zz+kk>0 && zz+kk<=n3
                               if dense(xx+ii,yy+jj,zz+kk)>15
                                 cn=cn+1;
                               end
                           end
                      end
                  end
             end
             if cn>5
                 dense_1(xx,yy,zz)=20;
             end
        end
    end
    end
end
            
%The processing factors could not run into the sphere of the center for
%pre-rRNA producing factory place
for x=1:n1
    for y=1:n2
        for z=1:n3
            if (x-(n1+1)/2)^2+(y-(n2+1)/2)^2+(z-(n3+1)/2)^2<=7
                position(x,y,z)=100;
            end
        end
    end
end

%calculation of the density of processing factors in the path from each
%processing factor producing place to the binding processing factor
for n=1:s(1)
    fac=randomFactor(n,:);
    inn=inner(fac(5),:);
    inn=round(inn);
    %xd=abs(inn(1)-fac(2));yd=abs(inn(2)-fac(2));zd=abs(inn(3)-fac(3));
    xd=inn(1)-fac(2);yd=inn(2)-fac(3);zd=inn(3)-fac(4);
    [mx,idx]=max([abs(xd),abs(yd),abs(zd)]);
    if xd<0
        x1=inn(1);x2=fac(2);
    else
        x1=fac(2);x2=inn(1);
    end
    if yd<0
        y1=inn(2);y2=fac(3);
    else
        y1=fac(3);y2=inn(2);
    end
    if zd<0
        z1=inn(3);z2=fac(4);
    else
        z1=fac(4);z2=inn(3);
    end
    
    count1=0;
    if idx==1
        for k=x1:x2
            density=sum(sum(sum((position(k,y1:y2,z1:z2)))));
            count1=count1+log(density+1);
        end
    end
    if idx==2
        for k=y1:y2
            density=sum(sum(sum(position(x1:x2,k,z1:z2))));
            count1=count1+log(density+1);
        end
    end
    if idx==3
        for k=z1:z2
            density=sum(sum(sum(position(x1:x2,y1:y2,k))));
            count1=count1+log(density+1);
        end
    end
    
    
   point2=fac(2:4);
   x2=point2(1);y2=point2(2);z2=point2(3);
   l1=abs(point2(1)-1);
   l2=abs(point2(1)-n1);
   l3=abs(point2(2)-1);
   l4=abs(point2(2)-n2);
   l5=abs(point2(3)-1);
   l6=abs(point2(3)-n3);
   [mn,idx]=min([l1,l2,l3,l4,l5,l6]);
   count2=0;
   if idx==1
       ys=y2-range;zs=z2-range;ye=y2+range;ze=z2+range;
       if ys<1 ys=1;end
       if ys>n2 ys=n2;end
       if zs<1 zs=1;end
       if zs>n3 zs=n3; end
       if ye<1 ye=1;end
       if ye>n2 ye=n2;end
       if ze<1 ze=1; end
       if ze>n3 ze=n3;end
       for kk=1:x2
            density=sum(sum(sum(position(kk,ys:ye,zs:ze))));
            count2=count2+log(density+1);
       end
   end
   if idx==2
       ys=y2-range;zs=z2-range;ye=y2+range;ze=z2+range;
       if ys<1 ys=1;end
       if ys>n2 ys=n2;end
       if zs<1 zs=1;end
       if zs>n3 zs=n3; end
       if ye<1 ye=1;end
       if ye>n2 ye=n2;end
       if ze<1 ze=1; end
       if ze>n3 ze=n3;end
       for kk=x2:n1
            density=sum(sum(sum(position(kk,ys:ye,zs:ze))));
            count2=count2+log(density+1);
       end
   end
   
   if idx==3
       xs=x2-range;zs=z2-range;xe=x2+range;ze=z2+range;
       if xs<1 xs=1;end
       if xs>n1 xs=n1;end
       if zs<1 zs=1;end
       if zs>n3 zs=n3; end
       if xe<1 xe=1;end
       if xe>n1 xe=n1;end
       if ze<1 ze=1; end
       if ze>n3 ze=n3;end
       for kk=1:y2
            density=sum(sum(sum(position(xs:xe,kk,zs:ze))));
            count2=count2+log(density+1);
       end
   end
    
    if idx==4
       xs=x2-range;zs=z2-range;xe=x2+range;ze=z2+range;
       if xs<1 xs=1;end
       if xs>n1 xs=n1;end
       if zs<1 zs=1;end
       if zs>n3 zs=n3; end
       if xe<1 xe=1;end
       if xe>n1 xe=n1;end
       if ze<1 ze=1; end
       if ze>n3 ze=n3;end
       for kk=y2:n2
            density=sum(sum(sum(position(xs:xe,kk,zs:ze))));
            count2=count2+log(density+1);
       end
    end
         
    if idx==5
       xs=x2-range;ys=y2-range;xe=x2+range;ye=y2+range;
       if xs<1 xs=1;end
       if xs>n1 xs=n1;end
       if ys<1 ys=1;end
       if ys>n3 ys=n3; end
       if xe<1 xe=1;end
       if xe>n1 xe=n1;end
       if ye<1 ye=1; end
       if ye>n3 ye=n3;end
       for kk=1:z2
            density=sum(sum(sum(position(xs:xe,ys:ye,kk))));
            count2=count2+log(density+1);
       end
    end
   
   if idx==6
       xs=x2-range;ys=y2-range;xe=x2+range;ye=y2+range;
       if xs<1 xs=1;end
       if xs>n1 xs=n1;end
       if ys<1 ys=1;end
       if ys>n3 ys=n3; end
       if xe<1 xe=1;end
       if xe>n1 xe=n1;end
       if ye<1 ye=1; end
       if ye>n3 ye=n3;end
       for kk=z2:n3
            density=sum(sum(sum(position(xs:xe,ys:ye,kk))));
            count2=count2+log(density+1);
       end
   end
    
    point2=fac(2:4);
    count3=dense(point2(1),point2(2),point2(3));
    eff=count3/(count1+count2);
    randomFactor(n,7)=eff;
    randomFactor(n,8)=count3;
    randomFactor(n,9)=count1;
    randomFactor(n,10)=count2;
end

randomFactor_1=randomFactor;
end
