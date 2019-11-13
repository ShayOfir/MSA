%input vector of values YY of order of magnitude -2 to 10 
clear sigmag fg xg gg calibYY ug vg gg4 dgg4 ug4 vg4

nn=length(YY);mm1g=min(YY);MM1g=max(YY);
sigmag=[.125 .25 .50 1]*((MM1g-mm1g)/2);
for j=1:4
    sigg=sigmag(j);
    for i=1:1000
       xg(i,j)=mm1g+(MM1g-mm1g)*i/1000;
       fg(i,j)=sum(exp(-(YY-xg(i,j)).^2/(2*sigg^2)))/sigg;
    end
end
aver=mean(YY);
gg=mean(fg')';
gg4=gg(xg(:,4)<aver);
dgg4=diff(gg4);
[ug4 vg4]=min(dgg4);
if ug4<0 & vg4>0
[ug vg]=max(gg4(1:vg4));
% i=2;while gg(i)>gg(i-1),i=i+1;end,vg=i-1;
if vg==vg4 | vg==1
    factor1=1;disp('no calibration')
end
if vg>1 & vg<vg4
    mode1=mm1g+(MM1g-mm1g)*vg/1000;
    factor1=1/(1-mode1/aver);
end
else
    factor1=1;disp('no calibration')
end
%if mode1 >= aver,factor=1;disp('no calibration'),mode1=aver;end
%if vg==length(gg(xg(:,4)<aver)),factor=1;disp('no calibration'),end
%if vg==1,factor=1;disp('no calibration'),end

%figure(1),plot(mm+(MM-mm)*(1:1000)'/1000,[g f]),grid,zoom
%disp([aver mode1 aver*length(YY) factor])
calibYY=aver+factor1*(YY-aver);
%figure(2),plot(calibYY),grid,zoom
