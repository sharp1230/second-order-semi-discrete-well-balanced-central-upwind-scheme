clear
clc

%% initial constant
N = 1201;
p = -1;
q = 1;
dx = (q-p)/(N-1);
g = 9.81;
t = 0;
theta = 1.3;
x = p : dx : q ;

v = Vinitial_2(x,N);        % Sets initial flow velocity
for i = 1 : N
    b(i) = Binitial_2(x(i));        % Sets initial flow depth
end
for i = 1 : N
    h(i) = Hinitial_2(x(i),b(i),0.00001);        % Sets initial flow depth
end

U1(1:N+2) = zeros;   % Defines array size of future grid points
U2(1:N+2) = zeros;

h = [h(1) h h(N)];
v = [v(1) v v(N)];
b = [b(1) b b(N)];
w = h+b;
hu = h.*v;
w0 = w(1);
wN1 = w(N);
hu0 = hu(1);
huN1 = hu(N);
dt = 0.25*dx;

%% main
while t<0.7
    w_old = w;
    hu_old = hu;
    DeltaX = dx;
    for p = 1:3
        wbar = w;
        hubar =hu;
        
        %% minmod of w_x
        wx(1) = minmod ( theta*(w(1)-w(1))/dx ,  theta*(w(2)-w(1))/dx , (w(2)-w(1))/2/dx );
        for j = 2:N-1
            wx(j) = minmod ( theta*(w(j)-w(j-1))/dx ,  theta*(w(j+1)-w(j))/dx , (w(j+1)-w(j-1))/2/dx );
        end
        wx(N) = minmod ( theta*(w(N)-w(N-1))/dx ,  theta*(w(N)-w(N))/dx , (w(N)-w(N-1))/2/dx );
        
        %% minmod of hu_x
        hux(1) = minmod ( theta*(hu(1)-hu(1))/dx ,  theta*(hu(2)-hu(1))/dx , (hu(2)-hu(1))/2/dx );
        for j = 2:N-1
            hux(j) = minmod ( theta*(hu(j)-hu(j-1))/dx ,  theta*(hu(j+1)-hu(j))/dx , (hu(j+1)-hu(j-1))/2/dx );
        end
        hux(N) = minmod ( theta*(hu(N)-hu(N-1))/dx ,  theta*(hu(N)-hu(N))/dx , (hu(N)-hu(N-1))/2/dx );
        
        %% w_l, w_r, hu_l, hu_r, h
        for j = 1:N
            w_r(j)  = w(j) - dx/2 * wx(j);
            hu_r(j) = hu(j) - dx/2 * hux(j);
        end
        w_r(N+1) = w(N);
        hu_r(N+1) = hu_r(N);
        
        for j = 2:N+1
            w_l(j)  = w(j-1) + dx/2 * wx(j-1);
            hu_l(j) = hu(j-1) + dx/2 * hux(j-1);
        end
        w_l(1) = w_r(2);
        hu_l(1) = hu_r(2);
        
        for j = 1:N
            h(j) = w(j) - (b(j+1)+b(j))/2;
        end
        
        %% a_plus, a_minus        
        for j = 1:N+1
            a_plus(j) = max( max( hu_r(j)/(w_r(j)-b(j)) + sqrt(w_r(j)-b(j)) , hu_l(j)/(w_l(j)-b(j)) + sqrt(w_l(j)-b(j)) ) , 0 );
            a_minus(j) = min( min( hu_r(j)/(w_r(j)-b(j)) - sqrt(w_r(j)-b(j)) , hu_l(j)/(w_l(j)-b(j)) - sqrt(w_l(j)-b(j)) ) , 0 );
        end
        
%         %% Special one
%         Lamda1 = 1/(4*max( max(abs(a_plus)),max(abs(a_minus)) ));
%         dt = Lamda1*dx;
        
        %% initial F
        for j = 1:N+1
            Fw_plus(j)     = hu_r(j);
            Fw_minus(j)  = hu_l(j);
            Fhu_plus(j)    = hu_r(j)^2/(w_r(j)-b(j)) + (w_r(j)-b(j))^2/2;
            Fhu_minus(j) = hu_l(j)^2/(w_l(j)-b(j)) + (w_l(j)-b(j))^2/2;
            
            wstar(j)   = (a_plus(j)*w_r(j) - a_minus(j)*w_l(j) - ( Fw_plus(j) - Fw_minus(j) ) ) / (a_plus(j) - a_minus(j));  
            hustar(j)  = (a_plus(j)*hu_r(j) - a_minus(j)*hu_l(j) - ( Fhu_plus(j) - Fhu_minus(j) ) ) / (a_plus(j) - a_minus(j)); 
        end
        
        %% initial d
        for j = 1:N+1
            dw(j) = minmod2( w_r(j)-wstar(j) , wstar(j) - w_l(j) ) / ( a_plus(j)-a_minus(j) );
            dhu(j) = minmod2( hu_r(j)-hustar(j) , hustar(j) - hu_l(j) ) / ( a_plus(j)-a_minus(j) );
        end
        
        %% compute Fw, Fhu
        for j = 1:N+1
            Fw(j)  = ( a_plus(j)*Fw_minus(j) - a_minus(j)*Fw_plus(j) ) / ( a_plus(j)-a_minus(j) ) + a_plus(j)*a_minus(j)*( (w_r(j)-w_l(j))/( a_plus(j)-a_minus(j))-dw(j) );
            Fhu(j) = ( a_plus(j)*Fhu_minus(j) - a_minus(j)*Fhu_plus(j) ) / ( a_plus(j)-a_minus(j) ) + a_plus(j)*a_minus(j)*( (hu_r(j)-hu_l(j))/( a_plus(j)-a_minus(j))-dhu(j) );
        end
        
        %% compute S2
        for j = 1:N
            S2(j) = - ((w_l(j+1)-b(j+1))+(w_r(j)-b(j)))/2*(b(j+1)-b(j)) / dx;
        end
        
        %% 3 state 3 order R-K method
        if p==1
            for j = 1:N
                k1w(j) = -(Fw(j+1) - Fw(j))/dx;
                k1hu(j) = -(Fhu(j+1) - Fhu(j))/dx + S2(j);
                w(j) = w(j) + dt*k1w(j);
                hu(j) = hu(j) + dt*k1hu(j);
            end
        elseif p==2
            for j = 1:N
                k2w(j) = -(Fw(j+1) - Fw(j))/dx;
                k2hu(j) = -(Fhu(j+1) - Fhu(j))/dx + S2(j);
                w(j) = w(j) + dt/4*(k1w(j)+k2w(j));
                hu(j) = hu(j) + dt/4*(k1hu(j)+k2hu(j));
            end
         elseif p==3
            for j = 1:N
                k3w(j) = -(Fw(j+1) - Fw(j))/dx;
                k3hu(j) = -(Fhu(j+1) - Fhu(j))/dx + S2(j);
            end
        end
    end
    for j = 1:N
        w(j) = w_old(j) + dt/6*(k1w(j)+k2w(j)+4*k3w(j));
        hu(j) = hu_old(j) + dt/6*(k1hu(j)+k2hu(j)+4*k3hu(j));
    end
    
    t = t + dt
    x = linspace(-1,1,N);

    plot(x, w(2:N+1),'b')  % Plots water height aproximations for each timestep.
    pause(0.01)
end


plot(w)










