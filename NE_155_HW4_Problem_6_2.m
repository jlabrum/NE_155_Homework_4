function [ x_GS,iterations_GS ] = NE_155_HW4_Problem_6_2( n,x_0 )
%Solves the equation Ax=b using the Gauss Seidel iteritive solving method.
%   Inputs for this function are n which is the length of x_0 and x_0 if
%the first x_k solution in the iteritive solving method. 
%   Outputs for this function are the final solution vector, x_GS, within an error
%greater than 10^-6 and the number of iterations, iterations_GS, it took.
format long
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

%Use the Guass Seidel Iteritive method to solve for x_GS.
%x_k1=-(D+L)^(-1)*U*x_k+(D+L)^(-1)*b, where L+D+U=A and L is strictly lower
%triangular, U is strictly upper triangular and D=diag(A)
D=diag(diag(A));
    
x=A\b;
x_k=x_0;


iterations_GS=0;
L=tril(A)-D;
U=triu(A)-D;
x_k=(-1)*inv(D+L)*U*x_0+inv(D+L)*b;
e=norm(abs(x_k-x_0));

while e > 10^-6
    x_k=(-1)*inv(D+L)*U*x_0+inv(D+L)*b;
    e=norm(abs(x_k-x_0));
    x_0=x_k;
    iterations_GS=iterations_GS+1;
   
end

x_GS=x_k;
end

