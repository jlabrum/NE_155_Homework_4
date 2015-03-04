function [ x_J,iterations_J] = NE_155_HW4_Problem_6_1(n,x_0)
%Solves the equation Ax=b using the Jacobi iteritive solving method.
%   Inputs for this function are n which is the length of x_0 and x_0 if
%the first x_k solution in the iteritive solving method. 
%   Outputs for this function are the final solution vector, x_J, within an error
%greater than 10^-6 and the number of iterations, iterations_J, it took.

%Build A and b based off of n
A=zeros(n);
A(1,1)=4;
A(1,2)=-1;
A(n,n-1)=-1;
A(n,n)=4;

for j=2:n-1
        A(j,j)=4;
        A(j,j+1)=-1;
        A(j,j-1)=-1;
end

b=100*ones(n,1);

%Use the Jacobi Iteritive method to solve for x_J.
%x_k1=D^(-1)*(D-A)*x_k+D^(-1)*b, where D is the diag(D)
D=diag(diag(A));
    
x=A\b;
x_k=x_0;

e=norm(abs(x-x_k));
iterations_J=0;
while e > 10^-6
    x_k=inv(D)*(D-A)*x_0+inv(D)*b;
    x_0=x_k;
    iterations_J=iterations_J+1;
    e=norm(abs(x-x_k));
end

x_J=x_k;




end

