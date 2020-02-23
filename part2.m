x= 3;           % x dimension
y=2;            % y dimension
%%%%% ratio for L/W is 3/2

dx = .1;
dy = .1;
nx = x/dx;
ny = y/dy;
%%create matrices for the F and G values. 
G = sparse(nx*ny,nx*ny);        %large matrix with mostly 0's
F = zeros(nx*ny,1);
% iterate through each point and give it the proper value. 
for i = 1:nx
    for j = 1:ny
        
        n = i + (j - 1) * nx; 
        nym = i + (j - 2) * nx; 
        nyp = i + j * nx; 
        nxm = i - 1 + (j - 1) * nx; 
        nxp = i + 1 + (j - 1) * nx; 
        
        if i == 1
            G(n, n) = 1; 
            F(n) = 1;
        elseif i == nx
            G(n, n) = 1; 
            F(n) = 1;
        elseif j == 1
            G(n, n) = 1; 
        elseif j == ny
            G(n, n) = 1; 
        else
            G(n, n) = -4; 
            G(n, nxm) = 1; 
            G(n, nxp) = 1; 
            G(n, nym) = 1; 
            G(n, nyp) = 1; 
        end
    end
end
%create the voltmap as well the matrix for V, density on X and Y. 
Voltmap2 = zeros(nx, ny); 
Ex=zeros (nx, ny);
Ey=zeros (nx, ny); 
V = G\F; 
%numerical way to answer voltage map 
for i = 1:nx
    for j = 1:ny
        n = i + (j - 1)*nx;
        Voltmap2(i, j) = V(n); 
    end
end
%iterate through the density matrices and find the corresponding density. 
 for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = Voltmap2(i+1, j) - Voltmap2(i, j); 
        elseif i == nx
            Ex(i, j) = Voltmap2(i, j) - Voltmap2(i-1, j);
        else
            Ex(i, j) = (Voltmap2(i+1, j) - Voltmap2(i-1, j))*0.5;
        end
        
        if j == 1
            Ey(i, j) = Voltmap2(i, j+1) - Voltmap2(i, j); 
        elseif j == ny
            Ey(i, j) = Voltmap2(i, j) - Voltmap2(i, j-1);
        else
            Ey(i, j) = (Voltmap2(i, j+1) - Voltmap2(i, j-1))*0.5;
        end
    end
 end
%plot the electric field for X and Y. 
figure(1)
surf(Ex)
title('Electric Field Ex')
xlabel('x')
ylabel('y')

figure(2)
surf(Ey)
title('Electric Field Ey')
xlabel('x')
ylabel('y')
%create the sigma matrix of ones. could have chosen 0, or 1. 
sigma = ones(nx,ny);
%iterate through the points and change where the bottle neck is. 
for i = 1:nx
    for j = 1:ny
        if j <= (ny/3) || j >= (ny*2/3)
            if i >= (nx/3) && i <= (nx*2/3)
                sigma(i,j) = 10^-12;
            end
            
        end
    end
end
%find the total electric field
TotalE= sqrt(Ex.^2+ Ey.^2);
%J is the total electric field where sigma is. 
J=sigma .* TotalE;
% plot sigma map
figure(3);
mesh(sigma);
xlabel('X');
ylabel('Y');
title('Sigma Map');
%plot current density map
figure(4)
mesh(J)
title('Current Density Map')
xlabel('x')
ylabel('y')
%%%%%%%Part 2 B%%%%%%%%
current =0;
OldCurrent=0;
TotalCurrent=0;
TotalE=0 ;
G = sparse(nx*ny,nx*ny);
F = zeros(nx*ny,1);
J=0;
TotalE=0;
sigma=0;
for meshsize = 1:5
    for i = 1:nx
        for j = 1:ny

            n = i + (j - 1) * nx; 
            nym = i + (j - 2) * nx; 
            nyp = i + j * nx; 
            nxm = i - 1 + (j - 1) * nx; 
            nxp = i + 1 + (j - 1) * nx; 

            if i == 1
                G(n, n) = 1; 
                F(n) = 1;
            elseif i == nx
                G(n, n) = 1; 
                F(n) = 1;
            elseif j == 1
                G(n, n) = 1; 
            elseif j == ny
                G(n, n) = 1; 
            else
                G(n, n) = -4; 
                G(n, nxm) = 1; 
                G(n, nxp) = 1; 
                G(n, nym) = 1; 
                G(n, nyp) = 1; 
            end
        end
    end

    Voltmap2 = zeros(nx, ny); 
    Ex=zeros (nx, ny);
    Ey=zeros (nx, ny); 
    V = G\F; 

    for i = 1:nx
        for j = 1:ny
            n = i + (j - 1)*nx;
            Voltmap2(i, j) = V(n); 
        end
    end
    
    sigma = ones(nx,ny);
        for i = 1:nx
        for j = 1:ny
            if j <= (ny/3) || j >= (ny*2/3)
                if i >= (nx/3) && i <= (nx*2/3)
                    sigma(i,j) = 10^-12;
                end

            end
        end
    end
     for i = 1:nx
        for j = 1:ny
            if i == 1
                Ex(i, j) = Voltmap2(i+1, j) - Voltmap2(i, j); 
            elseif i == nx
                Ex(i, j) = Voltmap2(i, j) - Voltmap2(i-1, j);
            else
                Ex(i, j) = (Voltmap2(i+1, j) - Voltmap2(i-1, j))*0.5;
            end

            if j == 1
                Ey(i, j) = Voltmap2(i, j+1) - Voltmap2(i, j); 
            elseif j == ny
                Ey(i, j) = Voltmap2(i, j) - Voltmap2(i, j-1);
            else
                Ey(i, j) = (Voltmap2(i, j+1) - Voltmap2(i, j-1))*0.5;
            end
        end
     end
       TotalE= sqrt(Ex.^2 + Ey.^2);
       J=sigma .* TotalE;
       OldCurrent=current;
       current =sum(J,'All');
       figure(5)
       plot([meshsize-1 meshsize], [OldCurrent current])
       hold on
       title ('Current vs Mesh Size')
end
% for meshsize = 1:ny        
%    TotalE= Ex(meshsize,:)+ Ey(meshsize,:);
%    J=sigma .* TotalE;
%    figure(5)
%    OldCurrent=TotalCurrent;
%    current =sum(J,2);
%    TotalCurrent=sum(current);
%    plot([meshsize-2 meshsize], [OldCurrent TotalCurrent],'o')
%    hold on
%    title ('Current vs Mesh Size')
% end
   xlabel('Mesh Size')
   ylabel('Current')

%%%%%%%%%%%%%%%%%%%%Part C%%%%%%%%%%%%%%
% change the size of our boxes inside dimensions and see effect on current
% density. 
% Sigma
sigma1 = ones(nx,ny);
sigma2 = ones(nx,ny);
sigma3 = ones(nx,ny);
sigma4 = ones(nx,ny);

for i = 1:nx %Changing the lengths and widths
        for j = 1:ny
            if j <= (ny/3) || j >= (ny*2/5)
                if i >= (nx/5) && i <= (nx*2/3)
                    sigma1(i,j) = 10^-12;
                end
            end
        if j <= (ny/4) || j >= (ny*3/4)
            if i >= (nx/2) && i <= (nx*2/3)
                sigma2(i,j) = 10^-2;
            end 
        end
        if j <= (ny/3) || j >= (ny*2/3)
            if i >= (nx/3) && i <= (nx*2/3)
                sigma3(i,j) = 10^-2;
            end 
        end 
        if j <= (ny/4) || j >= (ny/2)
            if i >= (nx/4) && i <= (nx/2)
                sigma4(i,j) = 10^-2;
            end 
        end
        end
end 

figure(6);
subplot(2,2,1)
mesh(sigma1);
xlabel('X');
ylabel('Y');
title('Sigma1 Map');
subplot(2,2,2)
mesh(sigma2);
xlabel('X');
ylabel('Y');
title('Sigma2 Map');
subplot(2,2,3)
mesh(sigma3);
xlabel('X');
ylabel('Y');
title('Sigma3 Map');
subplot(2,2,4)
mesh(sigma4);
xlabel('X');
ylabel('Y');
title('Sigma4 Map');
%%%%%%%%%%%%%%%Part D%%%%%%%%%%
sumJ=0;
oldJ=0;
for s= 1e-12:0.01:0.9
x= 3;           % x dimension
y=2;            % y dimension
%%%%% ratio for L/W is 3/2

dx = .1;
dy = .1;
nx = x/dx;
ny = y/dy;
%%create matrices for the F and G values. 
G = sparse(nx*ny,nx*ny);
F = zeros(nx*ny,1);

for i = 1:nx
    for j = 1:ny
        
        n = i + (j - 1) * nx; 
        nym = i + (j - 2) * nx; 
        nyp = i + j * nx; 
        nxm = i - 1 + (j - 1) * nx; 
        nxp = i + 1 + (j - 1) * nx; 
        
        if i == 1
            G(n, n) = 1; 
            F(n) = 1;
        elseif i == nx
            G(n, n) = 1; 
            F(n) = 1;
        elseif j == 1
            G(n, n) = 1; 
        elseif j == ny
            G(n, n) = 1; 
        else
            G(n, n) = -4; 
            G(n, nxm) = 1; 
            G(n, nxp) = 1; 
            G(n, nym) = 1; 
            G(n, nyp) = 1; 
        end
    end
end

Voltmap2 = zeros(nx, ny); 
Ex=zeros (nx, ny);
Ey=zeros (nx, ny); 
V = G\F; 

for i = 1:nx
    for j = 1:ny
        n = i + (j - 1)*nx;
        Voltmap2(i, j) = V(n); 
    end
end
 for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = Voltmap2(i+1, j) - Voltmap2(i, j); 
        elseif i == nx
            Ex(i, j) = Voltmap2(i, j) - Voltmap2(i-1, j);
        else
            Ex(i, j) = (Voltmap2(i+1, j) - Voltmap2(i-1, j))*0.5;
        end
        
        if j == 1
            Ey(i, j) = Voltmap2(i, j+1) - Voltmap2(i, j); 
        elseif j == ny
            Ey(i, j) = Voltmap2(i, j) - Voltmap2(i, j-1);
        else
            Ey(i, j) = (Voltmap2(i, j+1) - Voltmap2(i, j-1))*0.5;
        end
    end
end
Sigma = ones(nx,ny);

for i = 1:nx
    for j = 1:ny
        if j <= (ny/3) || j >= (ny*2/3)
            if i >= (nx/3) && i <= (nx*2/3)
                Sigma(i,j) = s;
            end
            
        end
    end
end
oldJ=sumJ;
TotalE= Ex + Ey;
J=Sigma .* TotalE;
sumJ=sum(J,'All');
figure(7);
hold on
plot([s s], [oldJ sumJ])
%mesh(Sigma);
xlabel('Sigma');
ylabel('Current');
title('Current vs Sigma');
end