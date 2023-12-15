x = [0; 1; 2; 0.5; 1.5; 0; 1; 2];
y = [0; 0; 0; 1; 1; 2; 2; 2];

C = [ [1; 2; 3; 5; 8; 7; 6; 4] [2; 3; 5; 8; 7; 6; 4; 1]];
%C = [1 5];
figure;
plot(x,y,'kx','LineWidth',2); grid on; hold on;
DT = delaunayTriangulation(x,y,C);
triplot(DT,'r-','LineWidth',1);

TF = isInterior(DT);