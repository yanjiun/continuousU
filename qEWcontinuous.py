import scipy.integrate
import numpy as np
import scipy.ndimage
import scipy.optimize

# 
# 08-08-2012 working on qEW simulations...
# make random field spline
# make force zero 
# advance u's one i at a time from flat initial configuration?
# 
# find zeros of velocity
#

class Subject:
    def __init__(self):
        self._observers = []

    def attach(self, observer):
        if not observer in self._observers:
            self._observers.append(observer)

    def detach(self, observer):
        try:
            self._observers.remove(observer)
        except ValueError:
            pass

    def notify(self, modifier=None):
        for observer in self._observers:
            if modifier != observer:
                observer.update(self)


class qEWContinuous(Subject):
    def __init__(self,m,N):
        Subject.__init__(self)
        self.randmat = 0.01*np.random.randn(4*N,N)
        self.us = np.zeros(N)
        self.m = m
        self.w = 0.0
        self.length = N
        self.VelocityArray = np.zeros(N)

    def findZero(self,x):
        # find next zero that will change the sign of the local Velocity
        # assume that current velocity is greater than zero    

        Felharm = self.Felharm(self.us[x],x)
        u_int = int(np.fix(self.us[x])+1)        

        # find next integer that makes the total velocity negative
        while Felharm + self.randmat[u_int,x] >= 0:
            u_int += 1
        # find zero in between current height and current_u_int
        print u_int, x

        #print self.localVelocity(self.us[x],x)
        #print self.localVelocity(u_int,x)
        unextzero = scipy.optimize.brentq(self.localVelocity,self.us[x],u_int,args=(x,))
        # replace current us[x] with the zero
        self.us[x] = unextzero        

    def findMetaStableString(self):
        """
        go through all indices to find the metastable string
        which order do we go?  does it matter?
        also need to treat periodic boundary conditions properly...
        will the string be found properly if we don't go in order?
        how to minimize the number of times this has to be done?
        """
        # first advance all indices smaller than 0?
        #for i in range(0,self.length):
        #    self.findZero(i)
        # now check if all localvelocities are zeros
        # and continue finding all the zeros until 
        self.calculateVelocityArray()
        while (self.VelocityArray>1e-16).any():
            indices = np.where(self.VelocityArray>1e-16)[0]
            for i in indices:
                self.findZero(i)
            self.calculateVelocityArray()

    def increaseW(self,deltaw):
        self.w += deltaw        

    def localVelocity(self,ux,x):
        laplacian = (self.us[x+1]-2*ux+self.us[x-1])
        quenched_noise = self.RFSpline(x,ux)

        return self.m**2*(self.w-ux) + laplacian + quenched_noise

    def Felharm(self,ux,x):
        laplacian = (self.us[x+1]-2*ux+self.us[x-1])
        
        return self.m**2*(self.w-ux) + laplacian

    def calculateVelocityArray(self):

        laplacian = scipy.ndimage.convolve1d(self.us,[1,-2,1],mode='wrap')
        quenched_noise = self.RFSplineArray()
        self.VelocityArray = self.m**2*(self.w-self.us) + laplacian +quenched_noise

    def RFSplineArray(self):
        """
        this returns an array of Random Fields, length of u array
        """

        u_lowers = np.fix(self.us)
        u_uppers = np.fix(self.us)+1
        lower_fields = self.randmat.transpose()[zip(*enumerate(u_lowers))]
        upper_fields = self.randmat.transpose()[zip(*enumerate(u_uppers))]        
        return lower_fields +(self.us-u_lowers)*(upper_fields-lower_fields)

    def RFSpline(self,x,u):
        """
        u is continuous, find closest two u integers
        this is for single x, mainly for finding roots
        """
        u_lower = int(np.fix(u)) 
        u_upper = int(np.fix(u))+1
        
        return self.randmat[u_lower,x] + (u-u_lower)*(self.randmat[u_upper,x]-self.randmat[u_lower,x])
