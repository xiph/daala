function [p32, p16, p8, c32, c16, c8] = process_bstats(stats)

X=sum(stats);

c32=X(1:28);
X=X(29:end);
d32=X(1:28);
X=X(29:end);

c16=X(1:144);
X=X(145:end);
d16=X(1:144);
X=X(145:end);

c8=c16-d16;
d8=zeros(16,144);
d8(:)=X;

p32 = d32./c32;
sum(H(p32).*c32)/sum(c32)

p16 = d16./(1e-100+c16);
sum(H(p16).*c16)/sum(c32)

p8 = d8./(1e-100+c8);
sum(sum(-p8.*log2(1e-100+p8).*c8))/sum(c32)
