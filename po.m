els=0:10:10;
azs=0:.5:360;
frs=[3,6];
myeps=1e-6;
k=2*pi/.3*permute(frs,[3,4,1,2]);
rcs=zeros(length(azs),length(frs));
v=load('v.txt');
f=load('f.txt');
v=cat(3,v(f(:,1),:),v(f(:,2),:),v(f(:,3),:));
l=circshift(v,[0,0,-1])-v;
c=(circshift(v,[0,0,-1])+v)/2;
a=sqrt(sum(cross(l(:,:,1),l(:,:,2),2).^2,2))/2;
n=cross(l(:,:,1),l(:,:,2))./a/2;
for ee=1:length(els)
  for aa=1:length(azs)
    inc=-[cosd(els(ee))*cosd(azs(aa)),cosd(els(ee))*sind(azs(aa)),sind(els(ee))];
    li=n*inc'<-myeps&n*inc'>-1+myeps;
    pe=n*inc'<=-1+myeps;
    rcs(aa,:)=20*log10(abs(permute(sum(-n(li,:)*inc'.*sum(sum(n(li,:)*[0,-inc(3),inc(2);inc(3),0,-inc(1);-inc(2),inc(1),0].*l(li,:,:),2).*exp(-2j*sum(c(li,:,:).*inc.*k,2)).*sinc(k.*sum(inc.*l(li,:,:),2)/pi),3)./(1-(n(li,:)*inc').^2)/2/sqrt(pi),1)+sum(-1j*k.*exp(-2j*c(pe,:,1)*inc'.*k).*a(pe)/sqrt(pi),1),[3,4,1,2])));
  end
  data=[azs',rcs];
  save([mfilename,'_',num2str(els(ee)),'.txt'],'data','-ascii');
  plot(azs,rcs);
  legend(num2str(frs'));
  grid on;
  xlabel('azimuth/degree');
  ylabel('RCS/dBsm');
  title(['elevation = ',num2str(els(ee)),' degree']);
  saveas(gcf,[mfilename,'_',num2str(els(ee)),'.png']);
end
