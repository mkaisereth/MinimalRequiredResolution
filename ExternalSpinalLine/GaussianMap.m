% https://math.stackexchange.com/questions/1301176/when-are-the-eigenvalues-of-the-second-fundamental-form-equal-to-the-principal-c
% https://en.wikipedia.org/wiki/Differential_geometry_of_surfaces
function [gc,dV,dL] = GaussianMap(x,y,z)
[xu,xv]     =   gradient(x);
[xuu,xuv]   =   gradient(xu);
[xvu,xvv]   =   gradient(xv);
[yu,yv]     =   gradient(y);
[yuu,yuv]   =   gradient(yu);
[yvu,yvv]   =   gradient(yv);
[zu,zv]     =   gradient(z);
[zuu,zuv]   =   gradient(zu);
[zvu,zvv]   =   gradient(zv);
dV = cell(size(z));
dL = cell(size(z));
for i=1:(size(z,1))
    for j=1:(size(z,2))
        Xu          =   [xu(i,j) yu(i,j) zu(i,j)];
        Xv          =   [xv(i,j) yv(i,j) zv(i,j)];
        Xuu         =   [xuu(i,j) yuu(i,j) zuu(i,j)];
        Xuv         =   [xuv(i,j) yuv(i,j) zuv(i,j)];
        Xvv         =   [xvv(i,j) yvv(i,j) zvv(i,j)];
        E           =   dot(Xu,Xu);
        F           =   dot(Xu,Xv);
        G           =   dot(Xv,Xv);
        m           =   cross(Xu,Xv);
        n           =   m/sqrt(sum(m.*m));
        L           =   dot(Xuu,n);
        M           =   dot(Xuv,n);
        N           =   dot(Xvv,n);
        % first fundamental matrix
        I = [E F; F G];
        % second fundamental matrix
        II = [L M; M N];
        % shape operator
        SO = II*inv(I);
        [dV1,dL1] = eig(SO);
        % eigenvectors of shape operator
        dV{i,j} = dV1;

        if ~isreal(dV1)
            disp('complex');
        end

        % eigenvalues of shape operator
        dL{i,j} = diag(dL1);
        % gaussian curvature is determinant of Shape operator
        gc1 = det(SO);
        gc(i,j)      =   ((L*N)-M^2)/((E*G)-F^2);
    end
end