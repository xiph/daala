function [rate,psnr] = bjontegaard(rd1,rd2,type)

    function [p] = pwl_fit(x,y)
        n=size(x,1)-1;
        p=[diff(y),y(1:n,:)];
    end

    function [p] = cubic_fit(x,y)
        s=min(x);
        e=max(x);
        p=polyfit((x.-s)./(e-s),y,3);
    end

    function [p] = interp_helper(x,y,k)
        n=size(x,1)-1;
        p=zeros(n,4);
        for i=1:n
            dx=x(i+1)-x(i);
            dy=y(i+1)-y(i);
            a=k(i)*dx-dy;
            b=-k(i+1)*dx+dy;
            p(i,4)=y(i);
            p(i,3)=a+dy;
            p(i,2)=b-2*a;
            p(i,1)=a-b;
        end
    end

    % Three-point first-order approximation
    function [k] = average_fit(x,y)
        n=size(x,1);
        dx=diff(x);
        dy=diff(y);
        d=dy./(2*dx);
        % one-sided average
        k=[d;0]+[0;d];
        % set end derivatives to zero
        %k(1)=0;
        %k(n)=0;
    end

    % Three-point second-order approximation
    function [k] = threepoint_fit(x,y)
        n=size(x,1);
        dx=diff(x);
        dy=diff(y);
        d=dy./dx;
        k=zeros(n,1);
        for i=2:n-1
            k(i)=(dx(i-1)*d(i)+dx(i)*d(i-1))/(dx(i)+dx(i-1));
        end
        % set end derivatives to zero
        k(1)=0;
        k(n)=0;
    end

    % Cubic Spline Interpolation
    function [k] = spline_fit(x,y)
        n=size(x,1);
        dx=diff(x);
        dy=diff(y);
        m=zeros(n,n);
        d=zeros(n,1);
        for i=1:n-1
            m(i,i+1)=1/dx(i);
            m(i+1,i)=1/dx(i);
            m(i+1,i+1)=2*m(i,i+1);
            m(i,i)=m(i,i)+2*m(i,i+1);
            d(i+1)=3*dy(i)/(dx(i)^2);
            d(i)=d(i)+3*dy(i)/(dx(i)^2);
        end
        k=m\d;
    end

    % Shape-preserving Cubic Hermite Interpolation
    function [k] = spchi_fit(x,y)
        n=size(x,1);
        dx=diff(x);
        dy=diff(y);
        d=dy./dx;

        % Assume all interior points have slope 0
        k=zeros(n,1);

        % Find index of non-local extrema interior points
        i=find(sign(d(1:n-2)).*sign(d(2:n-1))>0);

        % Compute the shape-preserving, non-centered three point average
        %  using Brodlie modification of Butland formula
        ds=dx(i)+dx(i+1);
        w1=(dx(i)+ds)./(3*ds);
        w2=(ds+dx(i+1))./(3*ds);
        dmax=max(abs(d(i)),abs(d(i+1)));
        dmin=min(abs(d(i)),abs(d(i+1)));
        k(i+1)=dmin./conj(w1.*(d(i)./dmax)+w2.*(d(i+1)./dmax));

        % Set the end points
        k(1)=((2*dx(1)+dx(2))*d(1)-dx(1)*d(2))/(dx(1)+dx(2));
        if sign(k(1))~=sign(d(1))
            k(1)=0;
        elseif (sign(d(1))~=sign(d(2)))&&(abs(k(1))>abs(3*d(1)))
            k(1)=3*d(1);
        end

        k(n)=((2*dx(n-1)+dx(n-2))*d(n-1)-dx(n-1)*d(n-2))/(dx(n-1)+dx(n-2));
        if sign(k(n))~=sign(d(n-1))
            k(n)=0;
        elseif (sign(d(n-1))~=sign(d(n-2)))&&(abs(k(n))>abs(3*d(n-1)))
            k(n)=3*d(n-1);
        end
    end

    function k = monotonize(x,y,k)
        n=size(x,1);
        dx=diff(x);
        dy=diff(y);
        d=dy./dx;
        a=k(1:n-1)./d;
        b=k(2:n)./d;
        t=[a<0;0]+[0;b<0];
        for i=1:n-1
            if (dy(i)==0)
                k(i)=0;
                k(i+1)=0;
                a(i+1)=0;
                b(i)=0;
                t=[a<0;0]+[0;b<0];
            elseif (t(i)>0)
                k(i)=0;
            elseif (a(i)*a(i)+b(i)*b(i)>9)
                tau=3/sqrt(a(i)*a(i)+b(i)*b(i));
                k(i)=tau*a(i)*d(i);
                k(i+1)=tau*b(i)*d(i);
            end
        end
    end

    function h = plot_poly(x,p)
        n=100;

        u=zeros(1,size(p,1)*n+1);
        v=zeros(1,size(p,1)*n+1);
        t=0:1/n:1;

        s=1;
        for i=1:size(p,1),
            u(s:(s+n))=t*(x(i+1,1)-x(i,1))+x(i,1);
            v(s:s+n)=polyval(p(i,:),t);
            s=s+n;
        end

        h = plot(u,v);
        hold all;
    end

    function val = int_poly(x,p,min_x,max_x)
        n=size(x,1);
        val=0;
        for i=1:n-1
            if abs(x(i+1,1)-x(i,1))>eps
                if x(i,1)>x(i+1,1)
                    x0=x(i+1,1);
                    x1=x(i,1);
                    t0=1;
                    t1=0;
                else
                    x0=x(i,1);
                    x1=x(i+1,1);
                    t0=0;
                    t1=1;
                end
                if x1>=min_x&&x0<=max_x
                    dx=x1-x0;
                    dt=t1-t0;
                    if min_x>x0
                        t0=t0+dt*(min_x-x0)/dx;
                    end
                    if max_x<x1
                        t1=t1+dt*(max_x-x1)/dx;
                    end
                    pi=polyint(p(i,:));
                    val=val+(polyval(pi,t1)-polyval(pi,t0))*(x1-x0);
                end
            end
        end
    end

    function val = avg_diff(x1,y1,x2,y2,type)
        min_x=max(min(x1),min(x2));
        max_x=min(max(x1),max(x2));

        switch type
            case 1
                p1=pwl_fit(x1,y1);
                p2=pwl_fit(x2,y2);
            case 2
                p1=cubic_fit(x1,y1);
                p2=cubic_fit(x2,y2);
                x1=[min(x1);max(x1)];
                x2=[min(x2);max(x2)];
            case 3
                p1=interp_helper(x1,y1,monotonize(x1,y1,average_fit(x1,y1)));
                p2=interp_helper(x2,y2,monotonize(x2,y2,average_fit(x2,y2)));
            case 4
                p1=interp_helper(x1,y1,monotonize(x1,y1,threepoint_fit(x1,y1)));
                p2=interp_helper(x2,y2,monotonize(x2,y2,threepoint_fit(x2,y2)));
            case 5
                p1=interp_helper(x1,y1,spline_fit(x1,y1));
                p2=interp_helper(x2,y2,spline_fit(x2,y2));
            case 6
                p1=interp_helper(x1,y1,spchi_fit(x1,y1));
                p2=interp_helper(x2,y2,spchi_fit(x2,y2));
        end

        int1=int_poly(x1,p1,min_x,max_x);
        int2=int_poly(x2,p2,min_x,max_x);

        val = (int2-int1)/(max_x-min_x);
    end

    rate1=log(rd1(:,1));
    psnr1=rd1(:,2);

    rate2=log(rd2(:,1));
    psnr2=rd2(:,2);

    psnr=avg_diff(rate1,psnr1,rate2,psnr2,type);
    rate=(exp(avg_diff(psnr1,rate1,psnr2,rate2,type))-1)*100;

    %p_pwl=pwl_fit(x1,y1);
    %p_cubic=cubic_fit(x1,y1);
    %p_average=interp_helper(x1,y1,monotonize(x1,y1,average_fit(x1,y1)));
    %p_threepoint=interp_helper(x1,y1,monotonize(x1,y1,threepoint_fit(x1,y1)));
    %p_spline=interp_helper(x1,y1,spline_fit(x1,y1));
    %p_spchi=interp_helper(x1,y1,spchi_fit(x1,y1));

    %plot_poly(x1,p_pwl);
    %plot_poly([min(x1);max(x1)],p_cubic);
    %plot_poly(x1,p_average);
    %plot_poly(x1,p_threepoint);
    %plot_poly(x1,p_spline);
    %plot_poly(x1,p_spchi);

    %int_poly(x1,p_pwl,min(x1),max(x1))
    %int_poly([min(x1);max(x1)],p_cubic,min(x1),max(x1))
    %int_poly(x1,p_average,min(x1),max(x1))
    %int_poly(x1,p_threepoint,min(x1),max(x1))
    %int_poly(x1,p_spline,min(x1),max(x1))
    %int_poly(x1,p_spchi,min(x1),max(x1))
    %plot(x1,y1,'o');
end
