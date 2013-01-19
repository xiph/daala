function [f] = sparsify(r)

    function [pgmi] = pgUpdate(Cn,Dn,E,F,SX,SY,m,i),
        ind=m(i,:)==1;
        Bni=Cn(ind,ind)\Dn(ind,i);
        B1i=SX(ind,ind)*Bni*SY(i,i)^-1;
        pgi=10*log10(F(i,i)/(F(i,i)-E(i,ind)*B1i));
        pgmi=zeros(1,64);
        for j=find(m(i,:)==1),
            m(i,j)=0;
            ind=m(i,:)==1;
            Bni=Cn(ind,ind)\Dn(ind,i);
            B1i=SX(ind,ind)*Bni*SY(i,i)^-1;
            pgmi(j)=pgi-10*log10(F(i,i)/(F(i,i)-E(i,ind)*B1i));
            m(i,j)=1;
        end
    end

    function [m,B1,B0] = fitModel(s,X,Y,m,d)
        Xm=mean(X);
        Xz=bsxfun(@minus,X,Xm);
        Ym=mean(Y);
        Yz=bsxfun(@minus,Y,Ym);

        C=Xz'*Xz;
        D=Xz'*Yz;
        E=Yz'*Xz;
        F=Yz'*Yz;

        SX=diag(1./diag(C).^(1/2));
        SY=diag(1./diag(F).^(1/2));

        Cn=SX*C*SX;
        Dn=SX*D*SY;

        pgm=zeros(16,64);

        % build the Pg masked look-up table
        for i=1:16,
            pgm(i,:)=pgUpdate(Cn,Dn,E,F,SX,SY,m,i);
        end

        % drop the top-d coefficients with the least impact on Pg
        for k=1:d,
            [pg,i]=min(pgm(m(:)==1));
            [i,j]=ind2sub(size(pgm),find(m(:)==1)(i));
            fprintf('Dropping (%i,%i) with pg=%g\n',i,j,pg)
            fflush(stdout);
            m(i,j)=0;
            pgm(i,:)=pgUpdate(Cn,Dn,E,F,SX,SY,m,i);
        end

        B1=zeros(64,16);
        B0=zeros(1,16);

        for i=1:16,
            ind=m(i,:)==1;
            Bni=Cn(ind,ind)\Dn(ind,i);
            B1(ind,i)=SX(ind,ind)*Bni*SY(i,i)^-1;
            B0(i)=Ym(i)-Xm(ind)*B1(ind,i);
        end

        fprintf('%s Blocks %i Pg=%g\n',s,size(X,1),sum(10*log10(diag(F)./diag(F-E*B1)))/16);
	fflush(stdout);
    end

    % reclassify based on weighted SATD
    function [c] = reclassify(X,Y,B1,B0)
        global OD_SCALE=[0.687397666275501251,0.691608410328626633,0.877750061452388763,0.874039031565189362];
        global scale=reshape(OD_SCALE'*OD_SCALE,1,16);

        fprintf('Reclassifying');
	fflush(stdout);
        for i=1:size(X,1),
            wSATD=zeros(10,1);
            if mod(i,10000)==0,
                fprintf('.');
                fflush(stdout);
            end
            for j=1:10,
                wSATD(j)=sum(abs(Y(i,:)-(B0(j,:)+X(i,:)*squeeze(B1(j,:,:)))).*scale);
            end
	    [v,c(i)]=min(wSATD);
	end
	fprintf('\n');
	fflush(stdout);
    end

    function [] = printMode(X,Y,i,B1,B0)
        Ym=mean(Y);
        Yz=bsxfun(@minus,Y,Ym);
	F=Yz'*Yz;
        E=Y-bsxfun(@plus,X*B1,B0);
        Em=mean(E);
        Ez=bsxfun(@minus,E,Em);
	G=Ez'*Ez;
        fprintf('Mode %i Blocks %i Pg=%g\n',i-1,size(X,1),sum(10*log10(diag(F)./diag(G))/16));
        fflush(stdout);
    end

    % 
    function [c,m,B1,B0] = kStep(c,X,Y,m,s,d)
        B1=zeros(10,64,16);
        B0=zeros(10,16);

        fprintf('Step %i (%i mults / block)\n',s,sum(m(:))/10-d);
	fflush(stdout);
        % fit the prediction model for 10 modes
        for i=1:10,
            ind=find(c(:)==i);
            [m(i,:,:),B1(i,:,:),B0(i,:)]=fitModel(sprintf('Mode %i',i-1),X(ind,:),Y(ind,:),squeeze(m(i,:,:)),d);
        end
	c=reclassify(X,Y,B1,B0);
	% print mode
        for i=1:10,
            ind=find(c(:)==i);
	    printMode(X(ind,:),Y(ind,:),i,squeeze(B1(i,:,:)),B0(i,:))
        end
    end

    c=r(:,3).+1;
    X=r(:,4:67);
    Y=r(:,68:83);

    m=ones(10,16,64);
    B1=zeros(10,64,16);
    B0=zeros(10,16);

    % initial 10 steps of full k-means
    for s=1:10,
        [c,m,B1,B0]=kStep(c,X,Y,m,s,0);
    end
    % drop from 1024 multiplies per mode to 64
    for s=11:70,
        [c,m,B1,B0]=kStep(c,X,Y,m,s,16);
    end
end
