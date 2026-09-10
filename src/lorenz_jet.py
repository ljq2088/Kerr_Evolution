"""Truncated two-variable Taylor algebra for local GHP operator composition.

Coefficients are derivatives divided by factorials in real r and theta.
No finite differencing is used: differential operators consume Taylor orders.
"""
import math
import numpy as np
from scipy.signal import convolve2d


class Jet:
    __array_priority__=1000

    def __init__(self, value=0., order=6, coefficients=None):
        self.order=order
        self.c=np.zeros((order+1,order+1),complex)
        if coefficients is None:
            self.c[0,0]=value
        else:
            self.c[:]=coefficients
        self.c[np.add.outer(np.arange(order+1),np.arange(order+1))>order]=0

    @classmethod
    def variable(cls,value,axis,order=6):
        out=cls(value,order)
        out.c[1,0] = int(axis==0)
        out.c[0,1] = int(axis==1)
        return out

    @property
    def value(self):
        return self.c[0,0]

    def lift(self,other):
        if isinstance(other,Jet):
            if other.order!=self.order:
                raise ValueError('Taylor orders differ')
            return other
        return Jet(other,self.order)

    def __add__(self,other):
        return Jet(order=self.order,coefficients=self.c+self.lift(other).c)
    __radd__=__add__

    def __neg__(self):
        return Jet(order=self.order,coefficients=-self.c)

    def __sub__(self,other):
        return self+-self.lift(other)

    def __rsub__(self,other):
        return self.lift(other)+-self

    def __mul__(self,other):
        rhs=self.lift(other)
        c=convolve2d(self.c,rhs.c)[:self.order+1,:self.order+1]
        return Jet(order=self.order,coefficients=c)
    __rmul__=__mul__

    def __pow__(self,power):
        if isinstance(power,int) and power>=0:
            result=Jet(1,self.order)
            for _ in range(power):
                result=result*self
            return result
        if self.value==0:
            raise ValueError('Noninteger or negative power at a zero expansion point')
        reduced=(self-self.value)/self.value
        term=Jet(1,self.order)
        result=term
        for n in range(1,self.order+1):
            term=term*reduced*(power-n+1)/n
            result=result+term
        return self.value**power*result

    def __truediv__(self,other):
        if not isinstance(other,Jet):
            return Jet(order=self.order,coefficients=self.c/other)
        return self*(other**-1)

    def __rtruediv__(self,other):
        return other*(self**-1)

    def derivative(self,axis):
        c=np.zeros_like(self.c)
        if axis==0:
            c[:-1,:]=self.c[1:,:]*np.arange(1,self.order+1)[:,None]
        else:
            c[:,:-1]=self.c[:,1:]*np.arange(1,self.order+1)[None,:]
        return Jet(order=self.order,coefficients=c)

    def exp(self):
        reduced=self-self.value
        term=Jet(1,self.order)
        result=term
        for n in range(1,self.order+1):
            term=term*reduced/n
            result=result+term
        return np.exp(self.value)*result

    def sin(self):
        return ((1j*self).exp()-(-1j*self).exp())/(2j)

    def cos(self):
        return ((1j*self).exp()+(-1j*self).exp())/2

    def conjugate(self):
        return Jet(order=self.order,coefficients=self.c.conjugate())

    def derivative_value(self,nr,nt):
        if nr+nt>self.order:
            raise ValueError('Derivative exceeds retained Taylor order')
        return self.c[nr,nt]*math.factorial(nr)*math.factorial(nt)
