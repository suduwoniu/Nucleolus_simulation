function[newxyz]=rotateEllipsoid(vert,azel,alpha)
alph=alpha*pi/180;
cosa=cos(alph);
sina=sin(alph);
vera=1-cosa;
theta=pi*azel(1)/180;
phi=pi*azel(2)/180;
u=[cos(phi)*cos(theta);cos(phi)*sin(theta);sin(phi)];
x=u(1);y=u(2);z=u(3);
rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';
newxyz=[vert(:,1),vert(:,2),vert(:,3)]*rot;
end
