n=100;
A=zeros(n);
A(1,1)=2;
A(1,2)=-1;
A(100,99)=-1;
A(100,100)=2;

for j=2:n-1
        A(j,j)=2;
        A(j,j+1)=-1;
        A(j,j-1)=-1;
end

b=[0:99]';

CN=norm(A,2)*norm(inv(A),2);


x1=inv(A)*b;
x2=A\b;
y=1:100;
hold on
plot(x1,y,'g','linewidth',5)
plot(x2,y,'r')

title('Solving a System with Inverse Calculation vs Built in Opertator','FontSize',15)
xlabel('x','FontSize',15)
ylabel('y','FontSize',15)
legend('Inverse Calculation','Built in Operator')