# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:53:37 2021

@author: Abhay Singh Rawat
Title:Gauss-Legendre-Quadrature method for integration

"""
import scipy.integrate as integral
import numpy.polynomial.legendre as leg


#legendre
def P(n,x):
    if n==0:
        return 1
    elif n==1:
        return x
    else:
        return (((2 * n)-1)*x * P(n-1,x)-(n-1)*P(n-2, x))/n

#question:F(X)=(X^2)/(1+X^3)
def func(x):
    f=(x**2)/(1+x**3)
    return f

#(1-x^2)*P'n(x)=-n*x*Pn(x)+n*Pn-1(x)
#dp_dx=lambda n,x:(-n*x*P(n,x)+n*P(n-1,x))/(1-x**2)

print("enter limits")
a=float(input("a\t"))
b=float(input("b\t"))
#no. of steps
n=int(input("enter n"))


"""finding roots of nth legendre polynomial and storing it in t"""
Q=[]
for i in range (n):
    Q.append(0)
Q.append(1)
t=leg.legroots(Q)
print(t)

#solution finding
sol=0
for i in range(n):
    u=((b-a)/2)*t[i]+(a+b)/2
    #print(u)
    #weight
    w=2*(1-t[i]**2)/(((n+1)**2)*(P(n+1,t[i]))**2)
    #print(w)
    sol1=w*func(u)
    sol=sol+sol1   


#solution Found
main_sol=(b-a)*sol/2
print("the value of integral of the function with n=",n,"is\n",main_sol)

man_sol=integral.quadrature(func, 0, 3,maxiter=n,miniter=1)
man_sol=round(man_sol[0],10)
print("manual sol",man_sol)