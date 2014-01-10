function [p8, p16, p32] = process_bstats(stats)

X=sum(stats);

L32 = 28;
d32=zeros(3,L32);
d32(:)=X(1:L32*3);
X=X(L32*3+1:end);

L16 = 64;
d16=zeros(8,L16);
d16(:)=X(1:L16*8);
X=X(L16*8+1:end);

%c16=X(1:144);
%X=X(145:end);
%d16=X(1:144);
%c8=c16-d16;
%X=X(145:end);

d8=zeros(16,144);
d8(:)=X;

c32 = sum(d32);
p32 = d32./(1e-100+c32);

c16 = sum(d16);
p16 = d16./(1e-100+c16);

%p16 = d16./(1e-100+c16);
%sum(H(p16).*c16)/sum(c32)

c8 = sum(d8);
p8 = d8./(1e-100+c8);
%sum(sum(-p8.*log2(1e-100+p8).*c8))/sum(c32)
