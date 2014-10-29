%c=cfl4(:,2:2:end);
%l=cfl4(:,1:2:end);
l=cfl4(:,1:15);
c=cfl4(:,16:end);
ff = find((sqrt(mean(l'.^2))>100) .* (sqrt(mean(c'.^2))>100));
length(ff)
l2 = l(ff,:);
c2 = c(ff,:);
l2 = l;
c2 = c;
c2 = c2.*sign(sum(c2'.*l2')');
corr=sum(c2'.*l2')./sqrt(sum(c2'.*c2').*sum(l2'.*l2'));
ff=find(corr > .8);
length(ff)
l3 = l2(ff,:);
c3 = c2(ff,:);
reg=sum(l3.*c3)./sum(l3.*l3);
printf("%d, ", round(128*reg/max(reg)));
printf("\n")

%at this point, we need to use the zigzag to convert to
%raster order for intra.c
