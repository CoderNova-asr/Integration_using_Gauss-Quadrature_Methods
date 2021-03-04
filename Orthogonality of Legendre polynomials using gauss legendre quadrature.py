# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 00:53:37 2021
@author:CoderNova-asr

"""

import numpy.polynomial.legendre as leg


#legendre
def P(n,x):
    if n==0:
        return 1
    elif n==1:
        return x
    else:
        return (((2 * n)-1)*x * P(n-1,x)-(n-1)*P(n-2, x))/n

#question
def func(n,m,x):
    f=P(n,x)*P(m,x)
    return f

#(1-x^2)*P'n(x)=-n*x*Pn(x)+n*Pn-1(x)
#dp_dx=lambda n,x:(-n*x*P(n,x)+n*P(n-1,x))/(1-x**2)

#print("enter limits")
a=-1
b=1

"""something special only for you"""
print("We have a function:\nf(x)=P(n,x)*P(m,x)")
N=int(input("Enter n for the above function"))
M=int(input("Enter m for the above function"))
print("Degree of f(x)=",N+M)
#no. of steps
n=int(input("enter N so that degree of f(x) is smaller than equal to 2N-1\t"))


"""finding roots of nth legendre polynomial and storing it in t"""
Q=[]
for i in range (n):
    Q.append(0)
Q.append(1)
t=leg.legroots(Q)
#print(t)


#solution finding
sol=0
for i in range(n):
    u=((b-a)/2)*t[i]+(a+b)/2
    #print(u)
    #weight
    w=2*(1-t[i]**2)/(((n+1)**2)*(P(n+1,t[i]))**2)
    #print(w)
    sol1=w*func(N,M,u)
    sol=sol+sol1   


#solution Found
main_sol=(b-a)*sol/2
print("the value of integral of the function with n=",n,"is\n",float(main_sol))

print("\n\n\n*******Hence,********")
print("for \tn=m: integral is zero\n n!=m: integral=2/(2n+1)")
