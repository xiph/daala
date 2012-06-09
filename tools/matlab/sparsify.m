function [f] = sparsify(r)

    function [B,rms] = beta(x,y)
        n=size(x,1);
        m=size(x,2);
        x=x-ones(size(x,1),1)*mean(x);
        y=y-ones(size(y,1),1)*mean(y);
        C=x'*x;
        D=x'*y;
        E=y'*y;
        Sx=diag(1./diag(C).^(1/2));
        Sy=diag(1./diag(E).^(1/2));
        B=C\D;
        e=y-x*B;
        e(1,:);
        y(1,:);
        std(y)./std(e);
        B=Sx\B*Sy;
        rms=sqrt(sum(diag((e'*e)))/(n*m));
    end

    c=r(:,1);
    x=[r(:,2:37),r(:,42:45),r(:,50:53),r(:,58:61)];
    y=[r(:,38:41),r(:,46:49),r(:,54:57),r(:,62:65)];

    f=[];

    for j=0:9,
        j
        ind=find(c(:)==j);
        [B,rms]=beta(x(ind,:),y(ind,1));
        F=[rms];

        [A,index]=sort(B);
        for i=2:48,
            [B,rms]=beta(x(ind,index(i:48)),y(ind,1));
            F=[F;rms];
        end
        f=[f,F];
    end
end