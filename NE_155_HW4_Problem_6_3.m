function [ x_SOR,iterations_SOR ] = NE_155_HW4_Problem_6_3( n,x_0,w )
%Solves the equation Ax=b using the SOR iteritive solving method.
%   Inputs for this function are n which is the length of x_0, x_0 is
%the first x_k solution in the iteritive solving method and w is a weight 
%added to speed p iterations,0<w<2. 
%   Outputs for this function are the final solution vector, x_SOR, within an error
%greater than 10^-6 and the number of iterations, iterations_SOR, it took.
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

%Use the SOR Iteritive method to solve for x_SOR.
%x_k1=-(D+w*L)^(-1)*[(1-w)D-wU]*x_k+(D+w*L)^(-1)*b, where L+D+U=A and L is 
%strictly lower triangular, U is strictly upper triangular and D=diag(A)and
%w is a weight to increase iterations
D=diag(diag(A));
    
x=A\b;
x_k=x_0;


iterations_SOR=0;
L=tril(A)-D;
U=triu(A)-D;
x_k=inv(D+w*L)*((1-w)*D-w*U)*x_0+inv(D+w*L)*w*b;
e=norm(abs(x_k-x_0));

while e > 10^-6
    x_k=inv(D+w*L)*((1-w)*D-w*U)*x_0+inv(D+w*L)*w*b;
    e=norm(abs(x_k-x_0));
    x_0=x_k;
    iterations_SOR=iterations_SOR+1;   
end

x_SOR=x_k;

end

