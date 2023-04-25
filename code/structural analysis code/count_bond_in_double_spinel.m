function [bond] = count_bond_in_double_spinel(filename)
fid=fopen(filename,'r+');

title=fscanf(fid,'%s',1);
a0=fscanf(fid,'%f',1);
coord=zeros(3,3);
coord(1,:)=fscanf(fid,'%f',3);
coord(2,:)=fscanf(fid,'%f',3);
coord(3,:)=fscanf(fid,'%f',3);
sMg=fscanf(fid,'%s',1);sAl=fscanf(fid,'%s',1); sFe=fscanf(fid,'%s',1);sO=fscanf(fid,'%s',1);
nMg=fscanf(fid,'%d',1);nAl=fscanf(fid,'%d',1); nFe=fscanf(fid,'%d',1);nO=fscanf(fid,'%d',1);
coordtype=fscanf(fid,'%s',1);
XYZ_Mg=zeros(nMg,3);XYZ_Al=zeros(nAl,3);XYZ_Fe=zeros(nFe,3);XYZ_O=zeros(nO,3);
for i=1:nMg
    XYZ_Mg(i,:)=fscanf(fid,'%f',3);
end
for i=1:nAl
    XYZ_Al(i,:)=fscanf(fid,'%f',3);
end
for i=1:nFe
    XYZ_Fe(i,:)=fscanf(fid,'%f',3);
end
for i=1:nO
    XYZ_O(i,:)=fscanf(fid,'%f',3);
end
a0=1;
R_Mg_O=zeros(nO,nMg);
r=zeros(27,1);
d=zeros(27,3);
PBD_OMg_type=zeros(nO,nMg);
cOMg=zeros(nO*nMg,3);
t=1;
for i=1:nO
    for j=1:nMg
        r(1)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(2)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(3)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(4)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(5)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(6)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(7)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        r(8)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(9)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(10)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(11)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3))^2);
        r(12)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(13)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(14)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        r(15)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2))^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        r(16)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(17)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(18)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        r(19)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1))^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        r(20)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(21)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(22)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(23)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        r(24)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)^2);
        r(25)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        r(26)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        r(27)=sqrt((XYZ_O(i,1)-XYZ_Mg(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Mg(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)^2);
        d(1,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(2,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(3,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(4,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(5,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(6,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(7,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        d(8,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(9,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(10,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(11,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3))];
        d(12,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(13,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(14,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        d(15,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        d(16,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(17,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(18,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        d(19,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        d(20,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(21,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(22,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(23,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        d(24,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3)+a0)];
        d(25,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)+a0),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        d(26,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)+a0),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        d(27,:)=[(XYZ_O(i,1)-XYZ_Mg(j,1)-a0),(XYZ_O(i,2)-XYZ_Mg(j,2)-a0),(XYZ_O(i,3)-XYZ_Mg(j,3)-a0)];
        [rr,k]=min(r);
        R_Mg_O(i,j)=rr;
        PBD_OMg_type(i,j)=k;
        cOMg(t,:)=d(k,:);
        t=t+1;
    end
end

b_R_Mg_O=zeros(nO*nMg,1);
t=1;
for i=1:nO
    for j=1:nMg
        b_R_Mg_O(t)=R_Mg_O(i,j);
        t=t+1;
    end
end
b_R_Mg_O=sort(b_R_Mg_O);

R_Fe_O=zeros(nO,nFe);PBD_OFe_type=zeros(nO,nFe);
cOFe=zeros(nO*nFe,3);
t=1;
for i=1:nO
    for j=1:nFe
        r(1)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(2)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(3)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(4)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(5)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(6)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(7)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(8)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(9)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(10)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(11)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3))^2);
        r(12)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(13)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(14)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(15)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2))^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(16)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(17)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(18)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(19)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1))^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(20)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(21)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(22)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(23)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(24)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(25)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(26)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(27)=sqrt((XYZ_O(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)^2);
        d(1,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(2,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(3,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(4,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(5,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(6,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(7,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        d(8,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(9,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(10,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(11,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3))];
        d(12,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(13,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(14,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        d(15,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        d(16,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(17,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(18,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        d(19,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        d(20,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(21,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(22,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(23,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        d(24,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3)+a0)];
        d(25,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)+a0),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        d(26,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)+a0),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        d(27,:)=[(XYZ_O(i,1)-XYZ_Fe(j,1)-a0),(XYZ_O(i,2)-XYZ_Fe(j,2)-a0),(XYZ_O(i,3)-XYZ_Fe(j,3)-a0)];
        %r=[r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26];
        [rr,k]=min(r);
        R_Fe_O(i,j)=rr;
        PBD_OFe_type(i,j)=k;
        cOFe(t,:)=d(k,:);
        t=t+1;
    end
end


b_R_Fe_O=zeros(nO*nFe,1);t=1;
for i=1:nO
    for j=1:nFe
        b_R_Fe_O(t)=R_Fe_O(i,j);
        t=t+1;
    end
end
b_R_Fe_O=sort(b_R_Fe_O);



R_Al_O=zeros(nO,nFe);PBD_OAl_type=zeros(nO,nAl);
cOAl=zeros(nO*nAl,3);
t=1;
for i=1:nO
    for j=1:nAl
        r(1)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(2)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(3)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(4)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(5)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(6)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(7)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        r(8)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(9)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(10)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(11)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3))^2);
        r(12)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(13)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(14)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        r(15)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2))^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        r(16)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(17)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(18)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        r(19)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1))^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        r(20)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(21)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(22)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(23)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        r(24)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)+a0)^2);
        r(25)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        r(26)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        r(27)=sqrt((XYZ_O(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_O(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_O(i,3)-XYZ_Al(j,3)-a0)^2);
        d(1,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(2,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(3,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(4,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(5,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(6,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(7,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        d(8,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(9,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(10,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(11,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3))];
        d(12,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(13,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(14,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        d(15,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        d(16,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(17,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(18,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        d(19,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        d(20,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(21,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(22,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(23,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        d(24,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3)+a0)];
        d(25,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)+a0),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        d(26,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)+a0),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        d(27,:)=[(XYZ_O(i,1)-XYZ_Al(j,1)-a0),(XYZ_O(i,2)-XYZ_Al(j,2)-a0),(XYZ_O(i,3)-XYZ_Al(j,3)-a0)];
        %r=[r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26];
        [rr,k]=min(r);
        R_Al_O(i,j)=rr;
        PBD_OAl_type(i,j)=k;
        cOAl(t,:)=d(k,:);
        t=t+1;
    end
end


b_R_Al_O=zeros(nO*nAl,1);t=1;
for i=1:nO
    for j=1:nAl
        b_R_Al_O(t)=R_Al_O(i,j);
        t=t+1;
    end
end
b_R_Al_O=sort(b_R_Al_O);





R_Fe_Mg=zeros(nMg,nFe);PBD_MgFe_type=zeros(nMg,nFe);
for i=1:nMg
    for j=1:nFe
        r(1)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(2)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(3)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(4)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(5)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(6)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(7)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(8)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(9)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(10)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(11)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3))^2);
        r(12)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(13)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(14)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(15)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2))^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(16)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(17)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(18)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(19)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1))^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(20)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(21)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(22)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(23)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(24)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)+a0)^2);
        r(25)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(26)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        r(27)=sqrt((XYZ_Mg(i,1)-XYZ_Fe(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Fe(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Fe(j,3)-a0)^2);
        [rr,k]=min(r);
        R_Fe_Mg(i,j)=rr;
        PBD_MgFe_type(i,j)=k;
    end
end

R_Al_Mg=zeros(nMg,nAl);PBD_MgAl_type=zeros(nMg,nAl);
for i=1:nMg
    for j=1:nAl
        r(1)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(2)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(3)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(4)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(5)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(6)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(7)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        r(8)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(9)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(10)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(11)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3))^2);
        r(12)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(13)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(14)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        r(15)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2))^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        r(16)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(17)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(18)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        r(19)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1))^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        r(20)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(21)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(22)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(23)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        r(24)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)+a0)^2);
        r(25)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        r(26)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        r(27)=sqrt((XYZ_Mg(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Mg(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Mg(i,3)-XYZ_Al(j,3)-a0)^2);
        [rr,k]=min(r);
        R_Al_Mg(i,j)=rr;
        PBD_MgAl_type(i,j)=k;
    end
end

R_Al_Fe=zeros(nFe,nAl);PBD_FeAl_type=zeros(nFe,nAl);
for i=1:nFe
    for j=1:nAl
        r(1)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(2)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(3)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(4)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(5)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(6)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(7)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        r(8)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(9)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(10)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(11)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3))^2);
        r(12)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(13)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(14)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        r(15)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2))^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        r(16)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(17)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(18)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        r(19)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1))^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        r(20)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(21)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(22)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(23)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        r(24)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)+a0)^2);
        r(25)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)+a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        r(26)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)+a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        r(27)=sqrt((XYZ_Fe(i,1)-XYZ_Al(j,1)-a0)^2+(XYZ_Fe(i,2)-XYZ_Al(j,2)-a0)^2+(XYZ_Fe(i,3)-XYZ_Al(j,3)-a0)^2);
        [rr,k]=min(r);
        R_Al_Fe(i,j)=rr;
        PBD_FeAl_type(i,j)=k;
    end
end

Mg_O_Fe_1=zeros(nO,6);
t=1;
% c1=zeros(1,3);c2=zeros(1,3);
for i=1:nO
    for j=1:nMg
        for k=1:nFe
            if R_Mg_O(i,j)<0.255 && R_Fe_O(i,k)<0.255
                Mg_O_Fe_1(t,1)=i;
                Mg_O_Fe_1(t,2)=j;
                Mg_O_Fe_1(t,3)=k;
                Mg_O_Fe_1(t,4)=R_Mg_O(i,j);
                d1=R_Mg_O(i,j);
                Mg_O_Fe_1(t,5)=R_Fe_O(i,k);
                d2=R_Fe_O(i,k);
                c1=cOMg((i-1)*nMg+j,:);
                c2=cOFe((i-1)*nFe+k,:);
                p=(c1*c2')/(d1*d2);
                Mg_O_Fe_1(t,6)=real(acos(p)*180/pi);
                t=t+1;
            end
        end
    end
end
           
Mg_O_Fe_2=zeros(nO,6);
t=1;
for i=1:nO
    for j=1:nMg
        for k=1:nFe
            if R_Mg_O(i,j)>0.402 && R_Mg_O(i,j)<0.46 && R_Fe_O(i,k)<0.255
                Mg_O_Fe_2(t,1)=i;
                Mg_O_Fe_2(t,2)=j;
                Mg_O_Fe_2(t,3)=k;
                Mg_O_Fe_2(t,4)=R_Mg_O(i,j);
                d1=R_Mg_O(i,j);
                Mg_O_Fe_2(t,5)=R_Fe_O(i,k);
                d2=R_Fe_O(i,k);
                c1=cOMg((i-1)*nMg+j,:);
                c2=cOFe((i-1)*nFe+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Mg_O_Fe_2(t,6)=real(acos(p)*180/pi);
                    t=t+1;
                end
            end
        end
    end
end
%Mg_O_Fe_2(t,:)=[];

Mg_O_Fe_3=zeros(nO,6);
t=1;
for i=1:nO
    for j=1:nMg
        for k=1:nFe
            if R_Mg_O(i,j)<0.255 && R_Fe_O(i,k)>0.40 && R_Fe_O(i,k)<0.46
                Mg_O_Fe_3(t,1)=i;
                Mg_O_Fe_3(t,2)=j;
                Mg_O_Fe_3(t,3)=k;
                Mg_O_Fe_3(t,4)=R_Mg_O(i,j);
                d1=R_Mg_O(i,j);
                Mg_O_Fe_3(t,5)=R_Fe_O(i,k);
                d2=R_Fe_O(i,k);
                c1=cOMg((i-1)*nMg+j,:);
                c2=cOFe((i-1)*nFe+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Mg_O_Fe_3(t,6)=real(acos(p)*180/pi);
                    t=t+1;
                end
            end
        end
    end
end

Mg_O_Al_1=zeros(nO,6);
t=1;
% c1=zeros(1,3);c2=zeros(1,3);
for i=1:nO
    for j=1:nMg
        for k=1:nAl
            if R_Mg_O(i,j)<0.255 && R_Al_O(i,k)<0.255
                Mg_O_Al_1(t,1)=i;
                Mg_O_Al_1(t,2)=j;
                Mg_O_Al_1(t,3)=k;
                Mg_O_Al_1(t,4)=R_Mg_O(i,j);
                d1=R_Mg_O(i,j);
                Mg_O_Al_1(t,5)=R_Al_O(i,k);
                d2=R_Al_O(i,k);
                c1=cOMg((i-1)*nMg+j,:);
                c2=cOAl((i-1)*nAl+k,:);
                p=(c1*c2')/(d1*d2);
                Mg_O_Al_1(t,6)=acos(p)*180/pi;
                t=t+1;
            end
        end
    end
end
           
Mg_O_Al_2=zeros(nO,6);
t=1;
for i=1:nO
    for j=1:nMg
        for k=1:nAl
            if R_Mg_O(i,j)>0.402 && R_Mg_O(i,j)<0.46 && R_Al_O(i,k)<0.255
                Mg_O_Al_2(t,1)=i;
                Mg_O_Al_2(t,2)=j;
                Mg_O_Al_2(t,3)=k;
                Mg_O_Al_2(t,4)=R_Mg_O(i,j);
                d1=R_Mg_O(i,j);
                Mg_O_Al_2(t,5)=R_Al_O(i,k);
                d2=R_Al_O(i,k);
                c1=cOMg((i-1)*nMg+j,:);
                c2=cOAl((i-1)*nAl+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Mg_O_Al_2(t,6)=real(acos(p)*180/pi);
                    t=t+1;
                end
            end
        end
    end
end
%Mg_O_Al_2(t,:)=[];

Mg_O_Al_3=zeros(nO,6);
t=1;
for i=1:nO
    for j=1:nMg
        for k=1:nAl
            if R_Mg_O(i,j)<0.255 && R_Al_O(i,k)>0.40 && R_Al_O(i,k)<0.46
                Mg_O_Al_3(t,1)=i;
                Mg_O_Al_3(t,2)=j;
                Mg_O_Al_3(t,3)=k;
                Mg_O_Al_3(t,4)=R_Mg_O(i,j);
                d1=R_Mg_O(i,j);
                Mg_O_Al_3(t,5)=R_Al_O(i,k);
                d2=R_Al_O(i,k);
                c1=cOMg((i-1)*nMg+j,:);
                c2=cOAl((i-1)*nAl+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Mg_O_Al_3(t,6)=real(acos(p)*180/pi);
                    t=t+1;
                end
            end
        end
    end
end

Fe_O_Al_1=zeros(nO,6);
t=1;
% c1=zeros(1,3);c2=zeros(1,3);
for i=1:nO
    for j=1:nFe
        for k=1:nAl
            if R_Fe_O(i,j)<0.255 && R_Al_O(i,k)<0.255
                Fe_O_Al_1(t,1)=i;
                Fe_O_Al_1(t,2)=j;
                Fe_O_Al_1(t,3)=k;
                Fe_O_Al_1(t,4)=R_Fe_O(i,j);
                d1=R_Fe_O(i,j);
                Fe_O_Al_1(t,5)=R_Al_O(i,k);
                d2=R_Al_O(i,k);
                c1=cOFe((i-1)*nFe+j,:);
                c2=cOAl((i-1)*nAl+k,:);
                p=(c1*c2')/(d1*d2);
                Fe_O_Al_1(t,6)=acos(p)*180/pi;
                t=t+1;
            end
        end
    end
end
           
Fe_O_Al_2=zeros(nO,6);
t=1;
for i=1:nO
    for j=1:nFe
        for k=1:nAl
            if R_Fe_O(i,j)>0.402 && R_Fe_O(i,j)<0.46 && R_Al_O(i,k)<0.255
                Fe_O_Al_2(t,1)=i;
                Fe_O_Al_2(t,2)=j;
                Fe_O_Al_2(t,3)=k;
                Fe_O_Al_2(t,4)=R_Fe_O(i,j);
                d1=R_Fe_O(i,j);
                Fe_O_Al_2(t,5)=R_Al_O(i,k);
                d2=R_Al_O(i,k);
                c1=cOFe((i-1)*nFe+j,:);
                c2=cOAl((i-1)*nAl+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Fe_O_Al_2(t,6)=real(acos(p)*180/pi);
                    t=t+1;
                end
            end
        end
    end
end
%Fe_O_Al_2(t,:)=[];

Fe_O_Al_3=zeros(nO,6);
t=1;
for i=1:nO
    for j=1:nFe
        for k=1:nAl
            if R_Fe_O(i,j)<0.255 && R_Al_O(i,k)>0.40 && R_Al_O(i,k)<0.46
                Fe_O_Al_3(t,1)=i;
                Fe_O_Al_3(t,2)=j;
                Fe_O_Al_3(t,3)=k;
                Fe_O_Al_3(t,4)=R_Fe_O(i,j);
                d1=R_Fe_O(i,j);
                Fe_O_Al_3(t,5)=R_Al_O(i,k);
                d2=R_Al_O(i,k);
                c1=cOFe((i-1)*nFe+j,:);
                c2=cOAl((i-1)*nAl+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Fe_O_Al_3(t,6)=real(acos(p)*180/pi);
                    t=t+1;
                end
            end
        end
    end
end

Mg_O_Mg_1=zeros(nO,6);
R_MgO_1=zeros(nO,1);R_MgO_2=zeros(nO,1);
t=1;
for i=1:nO
    for j=1:nMg-1
        for k=j+1:nMg
            if (R_Mg_O(i,j)<0.255 && R_Mg_O(i,k)>0.40&& R_Mg_O(i,k)<0.46)||(R_Mg_O(i,k)<0.255 && R_Mg_O(i,j)>0.40&& R_Mg_O(i,j)<0.46)
                Mg_O_Mg_1(t,1)=i;
                Mg_O_Mg_1(t,2)=j;
                Mg_O_Mg_1(t,3)=k;
                Mg_O_Mg_1(t,4)=R_Mg_O(i,j);
                d1=R_Mg_O(i,j);
                Mg_O_Mg_1(t,5)=R_Mg_O(i,k);
                d2=R_Mg_O(i,k);
                c1=cOMg((i-1)*nMg+j,:);
                c2=cOMg((i-1)*nMg+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Mg_O_Mg_1(t,6)=real(acos(p)*180/pi);
                    if d1>d2
                        R_MgO_1(t)=d2;
                        R_MgO_2(t)=d1;
                    elseif d2>d1
                        R_MgO_1(t)=d1;
                        R_MgO_2(t)=d2;
                    end
                    t=t+1;
                end
            end
        end
    end
end

Mg_O_Mg_2=zeros(nO/2,6);
t=1;
for i=1:nO
    for j=1:nMg-1
        for k=j+1:nMg
            if R_Mg_O(i,j)<0.255 && R_Mg_O(i,k)<0.255
                Mg_O_Mg_2(t,1)=i;
                Mg_O_Mg_2(t,2)=j;
                Mg_O_Mg_2(t,3)=k;
                Mg_O_Mg_2(t,4)=R_Mg_O(i,j);
                d1=R_Mg_O(i,j);
                Mg_O_Mg_2(t,5)=R_Mg_O(i,k);
                d2=R_Mg_O(i,k);
                c1=cOMg((i-1)*nMg+j,:);
                c2=cOMg((i-1)*nMg+k,:);
                p=(c1*c2')/(d1*d2);
                Mg_O_Mg_2(t,6)=real(acos(p)*180/pi);
                t=t+1;
            end
        end
    end
end

Fe_O_Fe_1=zeros(nO,6);
t=1;
for i=1:nO
    for j=1:nFe-1
        for k=j+1:nFe
            if R_Fe_O(i,j)<0.255 && R_Fe_O(i,k)<0.255
                Fe_O_Fe_1(t,1)=i;
                Fe_O_Fe_1(t,2)=j;
                Fe_O_Fe_1(t,3)=k;
                Fe_O_Fe_1(t,4)=R_Fe_O(i,j);
                d1=R_Fe_O(i,j);
                Fe_O_Fe_1(t,5)=R_Fe_O(i,k);
                d2=R_Fe_O(i,k);
                c1=cOFe((i-1)*nFe+j,:);
                c2=cOFe((i-1)*nFe+k,:);

                p=(c1*c2')/(d1*d2);
                Fe_O_Fe_1(t,6)=real(acos(p)*180/pi);
                t=t+1;
            end
        end
    end
end
             
Fe_O_Fe_2=zeros(nO,6);
R_FeO_1=zeros(nO,1);R_FeO_2=zeros(nO,1);
t=1;
for i=1:nO
    for j=1:nFe-1
        for k=j+1:nFe
            if (R_Fe_O(i,j)<0.255 && R_Fe_O(i,k)>0.40 && R_Fe_O(i,k)<0.46)||(R_Fe_O(i,k)<0.255 && R_Fe_O(i,j)>0.40 && R_Fe_O(i,j)<0.46)
                Fe_O_Fe_2(t,1)=i;
                Fe_O_Fe_2(t,2)=j;
                Fe_O_Fe_2(t,3)=k;
                Fe_O_Fe_2(t,4)=R_Fe_O(i,j);
                d1=R_Fe_O(i,j);
                Fe_O_Fe_2(t,5)=R_Fe_O(i,k);
                d2=R_Fe_O(i,k);
                c1=cOFe((i-1)*nFe+j,:);
                c2=cOFe((i-1)*nFe+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Fe_O_Fe_2(t,6)=real(acos(p)*180/pi);
                    if d1>d2
                        R_FeO_1(t)=d2;
                        R_FeO_2(t)=d1;
                    elseif d2>d1
                        R_FeO_1(t)=d1;
                        R_FeO_2(t)=d2;
                    end
                    t=t+1;
                end
            end
        end
    end
end

Al_O_Al_1=zeros(nO,6);
t=1;
for i=1:nO
    for j=1:nAl-1
        for k=j+1:nAl
            if R_Al_O(i,j)<0.255 && R_Al_O(i,k)<0.255
                Al_O_Al_1(t,1)=i;
                Al_O_Al_1(t,2)=j;
                Al_O_Al_1(t,3)=k;
                Al_O_Al_1(t,4)=R_Al_O(i,j);
                d1=R_Al_O(i,j);
                Al_O_Al_1(t,5)=R_Al_O(i,k);
                d2=R_Al_O(i,k);
                c1=cOAl((i-1)*nAl+j,:);
                c2=cOAl((i-1)*nAl+k,:);

                p=(c1*c2')/(d1*d2);
                Al_O_Al_1(t,6)=real(acos(p)*180/pi);
                t=t+1;
            end
        end
    end
end
             
Al_O_Al_2=zeros(nO,6);
R_AlO_1=zeros(nO,1);R_AlO_2=zeros(nO,1);
t=1;
for i=1:nO
    for j=1:nAl-1
        for k=j+1:nAl
            if (R_Al_O(i,j)<0.255 && R_Al_O(i,k)>0.40 && R_Al_O(i,k)<0.46)||(R_Al_O(i,k)<0.255 && R_Al_O(i,j)>0.40 && R_Al_O(i,j)<0.46)
                Al_O_Al_2(t,1)=i;
                Al_O_Al_2(t,2)=j;
                Al_O_Al_2(t,3)=k;
                Al_O_Al_2(t,4)=R_Al_O(i,j);
                d1=R_Al_O(i,j);
                Al_O_Al_2(t,5)=R_Al_O(i,k);
                d2=R_Al_O(i,k);
                c1=cOAl((i-1)*nAl+j,:);
                c2=cOAl((i-1)*nAl+k,:);
                p=(c1*c2')/(d1*d2);
                if p<0.2588
                    Al_O_Al_2(t,6)=real(acos(p)*180/pi);
                    if d1>d2
                        R_AlO_1(t)=d2;
                        R_AlO_2(t)=d1;
                    elseif d2>d1
                        R_AlO_1(t)=d1;
                        R_AlO_2(t)=d2;
                    end
                    t=t+1;
                end
            end
        end
    end
end  

Mg_O_Fe_1_new=sortrows(Mg_O_Fe_1,6);
Mg_O_Fe_2_new=sortrows(Mg_O_Fe_2,6);
Mg_O_Fe_3_new=sortrows(Mg_O_Fe_3,6);
Mg_O_Al_1_new=sortrows(Mg_O_Al_1,6);
Mg_O_Al_2_new=sortrows(Mg_O_Al_2,6);
Mg_O_Al_3_new=sortrows(Mg_O_Al_3,6);
Fe_O_Al_1_new=sortrows(Fe_O_Al_1,6);
Fe_O_Al_2_new=sortrows(Fe_O_Al_2,6);
Fe_O_Al_3_new=sortrows(Fe_O_Al_3,6);
Mg_O_Mg_1_new=sortrows(Mg_O_Mg_1,6);
Mg_O_Mg_2_new=sortrows(Mg_O_Mg_2,6);
Fe_O_Fe_1_new=sortrows(Fe_O_Fe_1,6);
Fe_O_Fe_2_new=sortrows(Fe_O_Fe_2,6);
Al_O_Al_1_new=sortrows(Al_O_Al_1,6);
Al_O_Al_2_new=sortrows(Al_O_Al_2,6);

for i=length(Mg_O_Fe_2_new(:,1)):-1:1
    tmp=0;
    for j=1:length(Mg_O_Fe_1_new(:,1))
        if (Mg_O_Fe_2_new(i,2)==Mg_O_Fe_1_new(j,2)) && (Mg_O_Fe_2_new(i,3)==Mg_O_Fe_1_new(j,3)) && (Mg_O_Fe_2_new(i,6)<85)
            tmp=tmp+1;
        end
    end
    if tmp>0
        Mg_O_Fe_2_new(i,:)=[];
    end
end


for i=length(Mg_O_Al_2_new(:,1)):-1:1
    tmp=0;
    for j=1:length(Mg_O_Al_1_new(:,1))
        if (Mg_O_Al_2_new(i,2)==Mg_O_Al_1_new(j,2)) && (Mg_O_Al_2_new(i,3)==Mg_O_Al_1_new(j,3)) && (Mg_O_Al_2_new(i,6)<85)
            tmp=tmp+1;
        end
    end
    if tmp>0
        Mg_O_Al_2_new(i,:)=[];
    end
end

for i=length(Fe_O_Al_2_new(:,1)):-1:1
    tmp=0;
    for j=1:length(Fe_O_Al_1_new(:,1))
        if (Fe_O_Al_2_new(i,2)==Fe_O_Al_1_new(j,2)) && (Fe_O_Al_2_new(i,3)==Fe_O_Al_1_new(j,3)) && (Fe_O_Al_2_new(i,6)<85)
            tmp=tmp+1;
        end
    end
    if tmp>0
        Fe_O_Al_2_new(i,:)=[];
    end
end

Mg_O_Fe_1 = count_bond_number2(Mg_O_Fe_1_new, 0.25, 0.25);
Mg_O_Fe_2 = count_bond_number2(Mg_O_Fe_1_new, 0.2165, 0.2500);
Mg_O_Fe_3 = count_bond_number2(Mg_O_Fe_1_new, 0.2500, 0.2165);
Mg_O_Fe_4 = count_bond_number2(Mg_O_Fe_2_new, 0.4146, 0.2165);
Mg_O_Fe_5 = count_bond_number2(Mg_O_Fe_2_new, 0.4330, 0.2500);
Mg_O_Fe_6 = count_bond_number2(Mg_O_Fe_2_new, 0.4146, 0.2500);            
Mg_O_Fe_7 = count_bond_number2(Mg_O_Fe_2_new, 0.4330, 0.2165);
Mg_O_Fe_8 = count_bond_number2(Mg_O_Fe_3_new, 0.2165, 0.4146);
Mg_O_Fe_9 = count_bond_number2(Mg_O_Fe_3_new, 0.2500, 0.4330);
Mg_O_Fe_10 = count_bond_number2(Mg_O_Fe_3_new, 0.2500, 0.4146);
Mg_O_Fe_11 = count_bond_number2(Mg_O_Fe_3_new, 0.2165, 0.4330);

Mg_O_Al_1 = count_bond_number2(Mg_O_Al_1_new, 0.25, 0.25);
Mg_O_Al_2 = count_bond_number2(Mg_O_Al_1_new, 0.2165, 0.2500);
Mg_O_Al_3 = count_bond_number2(Mg_O_Al_1_new, 0.2500, 0.2165);
Mg_O_Al_4 = count_bond_number2(Mg_O_Al_2_new, 0.4146, 0.2165);
Mg_O_Al_5 = count_bond_number2(Mg_O_Al_2_new, 0.4330, 0.2500);
Mg_O_Al_6 = count_bond_number2(Mg_O_Al_2_new, 0.4146, 0.2500);            
Mg_O_Al_7 = count_bond_number2(Mg_O_Al_2_new, 0.4330, 0.2165);
Mg_O_Al_8 = count_bond_number2(Mg_O_Al_3_new, 0.2165, 0.4146);
Mg_O_Al_9 = count_bond_number2(Mg_O_Al_3_new, 0.2500, 0.4330);
Mg_O_Al_10 = count_bond_number2(Mg_O_Al_3_new, 0.2500, 0.4146);
Mg_O_Al_11 = count_bond_number2(Mg_O_Al_3_new, 0.2165, 0.4330);

Fe_O_Al_1 = count_bond_number2(Fe_O_Al_1_new, 0.25, 0.25);
Fe_O_Al_2 = count_bond_number2(Fe_O_Al_1_new, 0.2165, 0.2500);
Fe_O_Al_3 = count_bond_number2(Fe_O_Al_1_new, 0.2500, 0.2165);
Fe_O_Al_4 = count_bond_number2(Fe_O_Al_2_new, 0.4146, 0.2165);
Fe_O_Al_5 = count_bond_number2(Fe_O_Al_2_new, 0.4330, 0.2500);
Fe_O_Al_6 = count_bond_number2(Fe_O_Al_2_new, 0.4146, 0.2500);            
Fe_O_Al_7 = count_bond_number2(Fe_O_Al_2_new, 0.4330, 0.2165);
Fe_O_Al_8 = count_bond_number2(Fe_O_Al_3_new, 0.2165, 0.4146);
Fe_O_Al_9 = count_bond_number2(Fe_O_Al_3_new, 0.2500, 0.4330);
Fe_O_Al_10 = count_bond_number2(Fe_O_Al_3_new, 0.2500, 0.4146);
Fe_O_Al_11 = count_bond_number2(Fe_O_Al_3_new, 0.2165, 0.4330);


Mg_O_Mg_1 = count_bond_number2(Mg_O_Mg_1_new, 0.4146, 0.2165) + count_bond_number2(Mg_O_Mg_1_new, 0.2165, 0.4146);
Mg_O_Mg_2 = count_bond_number2(Mg_O_Mg_1_new, 0.4146, 0.2500) + count_bond_number2(Mg_O_Mg_1_new, 0.2500, 0.4146);
Mg_O_Mg_3 = count_bond_number2(Mg_O_Mg_1_new, 0.2500, 0.4300) + count_bond_number2(Mg_O_Mg_1_new, 0.4300, 0.2500);
Mg_O_Mg_6 = count_bond_number2(Mg_O_Mg_1_new, 0.2165, 0.4330) + count_bond_number2(Mg_O_Mg_1_new, 0.4330, 0.2165);
Mg_O_Mg_4 = count_bond_number2(Mg_O_Mg_2_new, 0.2500, 0.2500);
Mg_O_Mg_5 = count_bond_number2(Mg_O_Mg_2_new, 0.2165, 0.2500) + count_bond_number2(Mg_O_Mg_2_new, 0.2500, 0.2165);

Fe_O_Fe_1 = count_bond_number2(Fe_O_Fe_1_new, 0.2500, 0.2500);
Fe_O_Fe_2 = count_bond_number2(Fe_O_Fe_1_new, 0.2165, 0.2500) + count_bond_number2(Fe_O_Fe_1_new, 0.2500, 0.2165);
Fe_O_Fe_3 = count_bond_number2(Fe_O_Fe_2_new, 0.2165, 0.4146) + count_bond_number2(Fe_O_Fe_2_new, 0.4146, 0.2165);
Fe_O_Fe_4 = count_bond_number2(Fe_O_Fe_2_new, 0.2500, 0.4330) + count_bond_number2(Fe_O_Fe_2_new, 0.4330, 0.2500);
Fe_O_Fe_5 = count_bond_number2(Fe_O_Fe_2_new, 0.2500, 0.4146) + count_bond_number2(Fe_O_Fe_2_new, 0.4146, 0.2500);
Fe_O_Fe_6 = count_bond_number2(Fe_O_Fe_2_new, 0.2165, 0.4330) + count_bond_number2(Fe_O_Fe_2_new, 0.4330, 0.2165);

Al_O_Al_1 = count_bond_number2(Al_O_Al_1_new, 0.2500, 0.2500);
Al_O_Al_2 = count_bond_number2(Al_O_Al_1_new, 0.2165, 0.2500) + count_bond_number2(Al_O_Al_1_new, 0.2500, 0.2165);
Al_O_Al_3 = count_bond_number2(Al_O_Al_2_new, 0.2165, 0.4146) + count_bond_number2(Al_O_Al_2_new, 0.4146, 0.2165);
Al_O_Al_4 = count_bond_number2(Al_O_Al_2_new, 0.2500, 0.4330) + count_bond_number2(Al_O_Al_2_new, 0.4330, 0.2500);
Al_O_Al_5 = count_bond_number2(Al_O_Al_2_new, 0.2500, 0.4146) + count_bond_number2(Al_O_Al_2_new, 0.4146, 0.2500);
Al_O_Al_6 = count_bond_number2(Al_O_Al_2_new, 0.2165, 0.4330) + count_bond_number2(Al_O_Al_2_new, 0.4330, 0.2165);
bond = [Mg_O_Fe_1, Mg_O_Fe_2, Mg_O_Fe_3, Mg_O_Fe_4, Mg_O_Fe_5, Mg_O_Fe_6, Mg_O_Fe_7, Mg_O_Fe_8, Mg_O_Fe_9, Mg_O_Fe_10, Mg_O_Fe_11, Mg_O_Al_1, Mg_O_Al_2, Mg_O_Al_3, Mg_O_Al_4, Mg_O_Al_5, Mg_O_Al_6, Mg_O_Al_7, Mg_O_Al_8, Mg_O_Al_9, Mg_O_Al_10, Mg_O_Al_11, Fe_O_Al_1, Fe_O_Al_2, Fe_O_Al_3, Fe_O_Al_4, Fe_O_Al_5, Fe_O_Al_6, Fe_O_Al_7, Fe_O_Al_8, Fe_O_Al_9, Fe_O_Al_10, Fe_O_Al_11, Mg_O_Mg_1, Mg_O_Mg_2, Mg_O_Mg_3, Mg_O_Mg_4, Mg_O_Mg_5, Mg_O_Mg_6,Fe_O_Fe_1, Fe_O_Fe_2, Fe_O_Fe_3, Fe_O_Fe_4, Fe_O_Fe_5, Fe_O_Fe_6, Al_O_Al_1, Al_O_Al_2, Al_O_Al_3, Al_O_Al_4, Al_O_Al_5, Al_O_Al_6];
fclose all;
end

