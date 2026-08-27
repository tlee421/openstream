function[x_full,y,T_full,T_center,q_bot_full,q_bot_total,q_top_full,q_top_total,x,M2]=ftoglass_ann_optim(M,N,W,H,k,lh,q_flux,w_free,htc_wake,htc_free,htc_amb,T_inf,T_amb)
    %inputs:
        %M = [-] number of x-nodes
        %N = [-] number of y-nodes
        %W = [m] half width of wall
        %H = [m] thickness of wall
        %k = [W/m-K] thermal conductivity of 304SS source: https://asm.matweb.com/search/specificmaterial.asp?bassnum=mq304a
        %q_flux = [W/m^2] heat flux into fluid
        %w_free = [m] total width of free stream
        %htc_wake = [W/m^2-K] average HTC for wake
        %htc_free = [W/m^2-K] average HTC for free stream
        %htc_amb = [W/m^2-K] average HTC for ambient natural convection
        %T_inf = [K] temperature of fluid (average film temperature for single phase flow)
        %T_amb = [K] temperature of ambient environment     
    %outputs:
        %x_full = [m] full x-position matrix
        %y = [m] full y-position matrix
        %T_full = [C] full temperature matrix of wall

    dx = W/(M-1); %[m] width of x-node
    dy = H/(N-1); %[m] width of y-node 
    x = ((1:M)-1)*dx; %[m] x-position of node i
    y = ((1:N)-1)*dy; %[m] y-position of node j 
    [~, M2] = min(abs(w_free/2 - x)); %index location of node closest to half width of free stream
    C =spalloc(M*N,M*N,5*M*N);
    % average htc between wake and free stream (weighted by quantity of node coverage)
        htc_avg = (x(M2)-w_free/2+dx/2)/dx*htc_wake+abs((x(M2)-w_free/2-dx/2)/dx*htc_free);

%% Calculation of Temperature Matrix    
    % bottom left node
        C(1,1)=-(htc_amb*dx/2+k*dy/(2*dx)+k*dx/(2*dy));
        C(1,2)=k*dy/(2*dx);
        C(1,M+1)=k*dx/(2*dy);
        d(1,1)= -htc_amb*dx/2*T_amb;        
    % bottom boundary    
        C = C+sparse(2:(M-1),2:(M-1),-(htc_amb*dx+k*dy/dx+k*dx/dy),M*N,M*N);
        C = C+sparse(2:(M-1),(2:(M-1))+1,k*dy/(2*dx),M*N,M*N);
        C = C+sparse(2:(M-1),M+(2:(M-1)),k*dx/dy,M*N,M*N);
        C = C+sparse(2:(M-1),(2:(M-1))-1,k*dy/(2*dx),M*N,M*N);
        d(2:(M-1),1) = -htc_amb*dx*T_amb;    
    % bottom right corner
        C(M,M)=-(htc_amb*dx/2+k*dx/(2*dy)+k*dy/(2*dx));
        C(M,M+M)=k*dx/(2*dy);
        C(M,M-1)=k*dy/(2*dx);
        d(M,1)=-htc_amb*dx/2*T_amb;        
    % internal nodes
        i = repmat(2:(M-1),1,N-2);
        j=sortrows(repmat(2:(M-1),1,N-2)')';
        C = C+sparse(M*(j-1)+i,M*(j-1)+i,-(2*k*dx/dy+2*k*dy/dx),M*N,M*N);
        C = C+sparse(M*(j-1)+i,M*(j-2)+i,k*dx/dy,M*N,M*N);
        C = C+sparse(M*(j-1)+i,M*(j-1)+i+1,k*dy/dx,M*N,M*N);
        C = C+sparse(M*(j-1)+i,M*(j)+i,k*dx/dy,M*N,M*N);
        C = C+sparse(M*(j-1)+i,M*(j-1)+i-1,k*dy/dx,M*N,M*N);
        d(M*(j-1)+i,1)=0;            
    % left side nodes
        j = 2:(N-1);
        C = C+sparse(M*(j-1)+1,M*(j-1)+1,-(k*dx/dy+k*dy/dx),M*N,M*N);
        C = C+sparse(M*(j-1)+1,M*(j-2)+1,k*dx/(2*dy),M*N,M*N);
        C = C+sparse(M*(j-1)+1,M*(j-1)+2,k*dy/dx,M*N,M*N);
        C = C+sparse(M*(j-1)+1,M*(j)+1,k*dx/(2*dy),M*N,M*N);
        d(M*(j-1)+1,1) = 0;        
    % right side nodes
        j = 2:(N-1);
        C = C+sparse(M*(j-1)+M,M*(j-1)+M,-(k*dx/dy+k*dy/dx),M*N,M*N);
        C = C+sparse(M*(j-1)+M,M*(j)+M,k*dx/(2*dy),M*N,M*N);
        C = C+sparse(M*(j-1)+M,M*(j-1)+M-1,k*dy/dx,M*N,M*N);
        C = C+sparse(M*(j-1)+M,M*(j-2)+M,k*dx/(2*dy),M*N,M*N);
        d(M*(j-1)+M,1) = 0;        
    % top left corner
        if M2>2
            C(M*(N-1)+1,M*(N-1)+1)=-(k*dx/(2*dy)+k*dy/(2*dx)+htc_free*dx/2);
            C(M*(N-1)+1,M*(N-2)+1)=k*dx/(2*dy);
            C(M*(N-1)+1,M*(N-1)+2)=k*dy/(2*dx);
            d(M*(N-1)+1,1)=-(htc_free*dx/2*T_inf+q_flux*dx/2);
        else
            C(M*(N-1)+1,M*(N-1)+1)=-(k*dx/(2*dy)+k*dy/(2*dx)+htc_wake*dx/2);
            C(M*(N-1)+1,M*(N-2)+1)=k*dx/(2*dy);
            C(M*(N-1)+1,M*(N-1)+2)=k*dy/(2*dx);
            d(M*(N-1)+1,1)=-(htc_wake*dx/2*T_inf+q_flux*dx/2);
        end        
     % top free stream nodes
        i=2:(M2-1);
        C = C+sparse(M*(N-1)+i,M*(N-1)+i,-(k*dx/dy+k*dy/dx+htc_free*dx),M*N,M*N);
        C = C+sparse(M*(N-1)+i,M*(N-2)+i,k*dx/dy,M*N,M*N);
        C = C+sparse(M*(N-1)+i,M*(N-1)+i+1,k*dy/(2*dx),M*N,M*N);
        C = C+sparse(M*(N-1)+i,M*(N-1)+i-1,k*dy/(2*dx),M*N,M*N);
        d(M*(N-1)+i,1) = -(htc_free*dx*T_inf+q_flux*dx);        
     % top wake nodes
        i=(M2+1):(M-1);
        C = C+sparse(M*(N-1)+i,M*(N-1)+i,-(k*dx/dy+k*dy/dx+htc_wake*dx),M*N,M*N);
        C = C+sparse(M*(N-1)+i,M*(N-2)+i,k*dx/dy,M*N,M*N);
        C = C+sparse(M*(N-1)+i,M*(N-1)+i+1,k*dy/(2*dx),M*N,M*N);
        C = C+sparse(M*(N-1)+i,M*(N-1)+i-1,k*dy/(2*dx),M*N,M*N);
        d(M*(N-1)+i,1) = -(htc_wake*dx*T_inf+q_flux*dx);  
     % wake-freestream interface node (only exists if M2 is interior)
        if M2 < M
            C(M*(N-1)+M2,M*(N-1)+M2)=-(k*dx/dy+k*dy/dx+htc_avg*dx);
            C(M*(N-1)+M2,M*(N-2)+M2)=k*dx/dy;
            C(M*(N-1)+M2,M*(N-1)+M2+1)=k*dy/(2*dx);
            C(M*(N-1)+M2,M*(N-1)+M2-1)=k*dy/(2*dx);
            d(M*(N-1)+M2,1)=-(htc_avg*dx*T_inf+q_flux*dx);
        end
      % top right corner
        if M2>(M-2)
            C(M*N,M*N)=-(k*dx/(2*dy)+k*dy/(2*dx)+htc_free*dx/2);
            C(M*N,M*(N-2)+M)=k*dx/(2*dy);
            C(M*N,M*(N-1)+M-1)=k*dy/(2*dx);
            d(M*N,1)=-(htc_free*dx/2*T_inf+q_flux*dx/2);
        else
            C(M*N,M*N)=-(k*dx/(2*dy)+k*dy/(2*dx)+htc_wake*dx/2);
            C(M*N,M*(N-2)+M)=k*dx/(2*dy);
            C(M*N,M*(N-1)+M-1)=k*dy/(2*dx);
            d(M*N,1)=-(htc_wake*dx/2*T_inf+q_flux*dx/2);
        end       
    % solution for temperature matrix
        S = C\d;
        T = reshape(S,M,N)-273.15; %[C] temperature matrix of wall
    % matrices for full plate (rather than half)
        x_left  = W - x(end:-1:1);     % gives 0 → W
        x_right = W + x;               % gives W → 2W
        x_full  = [x_left, x_right];   % final domain: 0 → 0.036
        T_right  = T(end:-1:1, :);      % mirror right half
        T_left = T;                   % original is the left half
        T_full  = [T_left; T_right];   % stack vertically in x direction
    % calculation for external centerline temperature
        T_center = T(M,1);
        
  %% Calculation of Heat Transfer
    % heat transferred through bottom boundary
        q_bot(1,1) = htc_amb*(T_amb-(T(1,1)+273.15));
        q_bot(M,1) = htc_amb*(T_amb-(T(M,1)+273.15));
        q_bot(2:(M-1),1)= htc_amb*(T_amb-(T(2:(M-1),1)+273.15));
        q_bot_full = abs([q_bot; q_bot(end:-1:1,:)]);
        q_bot_total = sum(q_bot_full)*lh;
     % heat transferred through top boundary (energy balance check)
        q_top(1,1) = htc_free*(T_inf-(T(1,N)+273.15));
        q_top(M,1) = htc_wake*(T_inf-(T(M,N)+273.15));
        q_top(2:M2,1) = htc_free*(T_inf-(T(2:M2,N)+273.15));
        q_top((M2+1):(M-1),1) = htc_wake*(T_inf-(T((M2+1):(M-1),N)+273.15));
        q_top_full = abs([q_top; q_top(end:-1:1,:)]);
        q_top_total = sum(q_top_full)*lh; 
end
        
 
        
    
    
    
    
    
    
    
    
    
    
    
    
    