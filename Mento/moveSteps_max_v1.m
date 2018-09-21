%This function is for choosing the numbers of processing factors
%(select_num) to move to the place for higher binding ability place in the
%3D space; 
%max_point: the most highest binding ability place
%ratio: the ratio of number of processing factors/the number of pre-rRNA
%producing place
function [rand_sort_res,inner]=moveSteps_max_v1(tem_track,inner,step,select_num,max_point,n1,n2,n3,ratio,rrna)
%tem_space: record each processing factor's position in a 3D space
tem_space=zeros(n1,n2,n3);
si=size(tem_track);
for n=1:si(1)
    x=tem_track(n,2);y=tem_track(n,3);z=tem_track(n,4);
    tem_space(x,y,z)=1;
end
s=size(tem_track);
%rand_sort: sort the processing factors from the higher binding ability to
%the lower binding ability
rand_sort=sortrows(tem_track,-7);
rand_select=rand_sort(1:max_point,:);
rand_move=rand_sort((max_point+1):s(1),:);
%s=size(rand_move);
%p=s(1);
%p1=randperm(p);
%select processing factors to move
p2=(s(1)-max_point-select_num):(s(1)-max_point);

%select for the movement direction of the select_num processing factors
%with the lowest binding ability
for n=1:length(p2)
    x=rand_move(p2(n),2);y=rand_move(p2(n),3);z=rand_move(p2(n),4);
    %random choose max point
    %mp=randperm(max_point);
    dist_max=zeros(max_point,1);
    for nn=1:max_point
        xm=rand_select(nn,2);ym=rand_select(nn,3);zm=rand_select(nn,4);
        dist_max(nn)=sqrt((x-xm)^2+(y-ym)^2+(z-zm)^2);
    end
    [m,idx]=min(dist_max);
    xm=rand_select(idx,2);ym=rand_select(idx,3);zm=rand_select(idx,4);
    d=sqrt((x-xm)^2+(y-ym)^2+(z-zm)^2);
    vx=(xm-x)/d;vy=(ym-y)/d;vz=(zm-z)/d;
    vs=round(rand*step);
    move=[round(x+vx*vs),round(y+vy*vs),round(z+vz*vs)];
    count=0;
    while (move(1)<1 || move(1)>n1 || move(2)<1 || move(2)>n2 || move(3)<1 || move(3)>n3|| tem_space(move(1),move(2),move(3))~=0 || sqrt((move(1)-n1/2)^2+(move(2)-n2/2)^2+(move(3)-n3/2)^2)<rrna)
        vs=round(rand*step);
        move=[round(x+vx*vs),round(y+vy*vs),round(z+vz*vs)];
        count=count+1;
        if count>100
            xt=rand*step*2-step; yt=rand*step*2-step;zt=rand*step*2-step;
            move=[x+round(xt),y+round(yt),z+round(zt)];
        end
    end
    rand_move(p2(n),2)=move(1);rand_move(p2(n),3)=move(2);rand_move(p2(n),4)=move(3);
    rand_sort((p2(n)+max_point),2)=move(1);
    rand_sort((p2(n)+max_point),3)=move(2);
    rand_sort((p2(n)+max_point),4)=move(3);
    tem_space(move(1),move(2),move(3))=1;
    tem_space(x,y,z)=0;
end
    inner(:,4)=0;
    nprocess=s(1);
    s2=size(inner);
    nrna=s2(1);
    for nn=1:nprocess
        dist0=99999999999;idx=0;
        for k=1:nrna
            if inner(k,4)<ratio
                dist1=sqrt((rand_sort(nn,2)-inner(k,1))^2+(rand_sort(nn,3)-inner(k,2))^2+(rand_sort(nn,4)-inner(k,3))^2);
                if dist1<dist0
                    dist0=dist1;
                    idx=k;
                end
            end
        end
        rand_sort(nn,5)=idx;
        rand_sort(nn,6)=dist0;
        inner(idx,4)=inner(idx,4)+1;
    end
    rand_sort_res=CalEffect_all_new_v4(rand_sort,inner,n1,n2,n3,3);
end

