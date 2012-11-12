import scipy.integrate
import numpy as np
import scipy.ndimage
import scipy.optimize
import scipy.sparse
from scipy.sparse.linalg import spsolve
from scipy.linalg import eig, solve

# 
# 08-21-2012 working on qEW simulations...
# make matrix for linear problem you want to solve
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
    def __init__(self,m,N,gamma,sigma):
        Subject.__init__(self)
        self.randmat = sigma*np.random.randn(4*N,N)
        self.us = np.zeros(N)
        self.m = m
        self.w = 0.0
        self.length = N
        self.VelocityArray = np.zeros(N)
        self.lhs = np.zeros((N,N)) 
        self.gamma = gamma
        self.current_cs = 0
        self.current_as = 0    

    def buildMatrix(self):
        """
        this will only be called once per simulation?
        do i build this directly as a sparse matrix or as a dense matrix then convert
        make a self.lhs and self.rhs? 
        """
        gammaRow = np.zeros(self.length) 
        gammaRow[0] = self.gamma
        gammaRow[-1] = self.gamma
        a = gammaRow
        for i in range(0,self.length-1):
            np.roll(gammaRow,1)
            a=np.vstack((a,gammaRow)) 
        c_s, a_s = self.RFSplineCoeffList()
        self.current_cs = c_s
        self.current_as = a_s
        diag_terms = -self.m**2 - 2*self.gamma + c_s
        np.fill_diagonal(a,diag_terms)

        # now make the dense matrix sparse
        a = scipy.sparse.csc_matrix(a) #
        # set diagonal of a

        self.lhs = a
        self.rhs = -self.m**2*self.w - a_s

    def doesAvalancheHappen(self):
        """
        check dFi/duj matrix (d^2 E/ dui duj)
        if any of the eigenvalues are positive
        then an avalanche happens

        """
        # problem with line below! it is not possible to compute all 
        # eigenvectors of sparse matrix... need to convert back to dense matrix
        # this is probably both time and memory consuming -- think of way around it
        denselhs = self.lhs.todense()
        evals, evects = eig(denselhs) 
        self.evals = evals
        self.evects = evects

        if (evals >= 0).any():
            return True
        else:
            return False

    def solveForUs(self):
        """
        this is just to solve the equation...
        """
        us_solution = spsolve(self.lhs,self.rhs) 
      
        return us_solution

    def updateAllMatrix(self):
        """
        this updates the matrix once something has moved
        !!! actually if we keep a list of what has been moved we can update only those diagonal terms!!!
        time the two different method to see which one is faster...
        """
        c_s, a_s = self.RFSplineCoeffList()
        diag_terms = -self.m**2 - 2*self.gamma + c_s
        self.lhs.setdiag(diag_terms)
        self.rhs = -self.m**2*self.w - a_s
        self.current_cs = c_s
        self.current_as = a_s

    def updatePartialMatrix(self,list_moved):
        
        old_diagonal = self.lhs.diagonal()
        for i in list_moved:
            c, a = self.RFSplineCoeff(self.us[i],i) 
            old_diagonal[i] = -self.m**2 - 2*self.gamma + c
            self.rhs[i] = -self.m**2*self.w - a 
            self.current_cs[i] = c
            self.current_as[i] = a
        self.lhs.setdiag(old_diagonal)

    def findFirstMinPush(self):
        
        u_next = np.fix(self.us)+1
        mwpushes = self.m**2*(u_next-self.us) + 2*self.gamma*(u_next-self.us) -self.current_cs*(u_next-self.us)
        #maskedpushes = np.ma.array(mwpushes,mask=(mwpushes>0))
        #print maskedpushes
        print "local forces:", mwpushes
        minIndex = mwpushes.argmin()
        mwmin = mwpushes[minIndex]
            
        return mwmin, minIndex

    def findMinPush(self):
        """
        this is to find the general minimum push given a configuration
        """
        mwpush = []
        for i in range(0, self.length):
            u_next = np.fix(self.us[i])+1
            deltaF = self.localVelocity(u_next,i)-self.localVelocity(self.us[i],i) 
            # the above line can be simplified
            mwpush.append(deltaF)
        mwpush = np.array(mwpush)
        # don't allow backwards motion for now...
        maskedpushes = np.ma.array(mwpush,mask = (mwpush>0))
        minIndex = maskedpushes.argmin()
 
        return mwpush[minIndex], minIndex
        
    def checkForSitesInAvalanche(self, u_next ,mwpush):
        """
        should I evolve all sites with F>0 forward or mwpush?
        """         
        old_noise = self.RFSplineArray(self.us)
        new_noise = self.RFSplineArray(u_next)
        
        deltaF = new_noise - old_noise - self.m**2*(u_next - self.us)-2*self.gamma*(u_next-self.us)
        #sitesToPush = (deltaF<=mwpush)
        pushedarray = ((deltaF+mwpush)>0)
        sitesToPush = np.where(pushedarray)[0]

        # returns array of indices in avalanche.  advance those forward?
        return sitesToPush, pushedarray

    def findDeltaW(self):
        return self.VelocityArray.min()/self.m**2.

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


    def increaseW(self,deltaw):
        self.w += deltaw        

    def localVelocity(self,ux,x):
        # implementing for periodic boundary conditions for now...
        laplacian = (self.us[np.mod(x+1,self.length)]-2*ux+self.us[np.mod(x-1,self.length)])
        quenched_noise = self.RFSpline(x,ux)

        return self.m**2*(self.w-ux) + self.gamma*laplacian + quenched_noise

    def Felharm(self,ux,x):
        laplacian = (self.us[x+1]-2*ux+self.us[x-1])
        
        return self.m**2*(self.w-ux) + self.gamma*laplacian    


    def calculateVelocityArray(self, u_array):

        quenched_noise = self.RFSplineArray(u_array)
        
        return self.calculateFharmArray(u_array)+quenched_noise

    def calculateFharmArray(self,u_array):
        laplacian = scipy.ndimage.convolve1d(u_array,[1,-2,1],mode='wrap')
    
        return self.m**2*(self.w - u_array) + self.gamma*laplacian

    def RFSplineArray(self, u_array):
        """
        this returns an array of Random Fields, length of u array
        """

        u_lowers = np.fix(u_array)
        u_uppers = np.fix(u_array)+1
        lower_fields = self.randmat.transpose()[zip(*enumerate(u_lowers))]
        upper_fields = self.randmat.transpose()[zip(*enumerate(u_uppers))]        
        return lower_fields +(u_array-u_lowers)*(upper_fields-lower_fields)

    def RFSpline(self,x,u):
        """
        u is continuous, find closest two u integers
        this is for single x, mainly for finding roots
        """
        u_lower = int(np.fix(u)) 
        u_upper = int(np.fix(u))+1

        return self.randmat[u_lower,x] + (u-u_lower)*(self.randmat[u_upper,x]-self.randmat[u_lower,x])

    def RFSplineCoeffList(self):
        
        u_lowers = np.fix(self.us)
        u_uppers = np.fix(self.us)+1
        
        lower_fields = self.randmat.transpose()[zip(*enumerate(u_lowers))]
        upper_fields = self.randmat.transpose()[zip(*enumerate(u_uppers))]        
        c_s = upper_fields-lower_fields
        a_s = lower_fields - u_lowers*(upper_fields-lower_fields) 

        return c_s, a_s

    def RFSplineCoeff(self,u,x):
        
        u_lower = np.floor(u)
        u_upper = np.floor(u)+1
        
        lower_field = self.randmat[u_lower,x]
        upper_field = self.randmat[u_upper,x]
        
        c = upper_field-lower_field
        a = lower_field - u_lower*(upper_field-lower_field)

        return c, a

def main(m,N,gamma,sigma):
    
    testModel = qEWContinuous(m,N,gamma,sigma)

    # set initial external force equal to the most negative a_s
    c_s, a_s = testModel.RFSplineCoeffList() 
    #initialForce = abs(min(a_s))
    testModel.buildMatrix()
    initialForce, minIndex = testModel.findFirstMinPush()
    print initialForce, minIndex
    testModel.increaseW(initialForce/testModel.m**2)
    #testModel.increaseW(initialForce/testModel.m**2)    

    testModel.updateAllMatrix()
    usolution = testModel.solveForUs()

    print "initial solution:", usolution
    # want to push u forward until we are happy 
    # print warning if solutions are not within [0,1]
    # NOTE: initial solution needs to be within range of where we expect... 
    # how do i write this so that I can trickily ensure this happens!?!?!?!    
    # or for now let the solutions be between (-inf,1]... 

    if (usolution>1).any() or (usolution<0).any():
        print "initial configuration not within range..." 
    testModel.us = usolution.copy()
    print "found first solution"
    print testModel.us 
    
    # evolve configuration equentially
    # 1. find minW to push by calculating velocity array
    # 2. push that one over
    # 3. find all other sets that would move given the push and the movement of that? or i guess i should find the next one that moves, push that... and then so on...?
    mwpush, minIndex = testModel.findFirstMinPush() 
    testModel.us[minIndex] = np.floor(testModel.us[minIndex]) + 1
    u_next = np.fix(testModel.us)+1
    u_upper = np.fix(testModel.us)+1
    u_lower = np.fix(testModel.us)
    testModel.updatePartialMatrix([minIndex])
    # solve for the configuration here
    usolution = testModel.solveForUs()
    
    testModel.us = usolution

    # check if there is an unstable direction
    # need to find parameters where this does go unstable...
    print testModel.doesAvalancheHappen()
       
 
    # increase w a little bit then propagate avalanche by using the most unstable direction, until another site crosses a random force threshold, or a local force is zero...
     

    return testModel

"""  
    # check if avalanche happens given this push, if not, advance another site too   
    # this part needs to be in a while loop till the simulation reaches the end
    while (testModel.us < 4*testModel.length).all():
        if testModel.doesAvalancheHappen():
            
            #avalanche occured, so find where avalanche occurs
            
            avalancheEnd = False
            listmoved, pushedarray = testModel.checkForSitesInAvalanche(u_next, mwpush)
            while not avalancheEnd:
                testModel.updatePartialMatrix(listmoved)
                # see if solution stopped in range 
                usolution = testModel.solveForUs() 
                u_upper += 1.*(pushedarray)
                u_lower +=1.*(pushedarray)
                # check if solution is within proper range
                # if not, then push over sites that should be pushed over
                if (usolution < u_lower).any():
                    print "solution behind range expected"
                checksolution = (usolution<=u_upper)
                if not checksolution.all():
                    # where check solution is false should be pushed over too
                    pushedarray = (checksolution*1)==0
                    listmoved = np.where(pushedarray)[0]
                else:
                    #avalanche ended!, update also testModel.us with solution 
                    avalancheEnd = True
                    testModel.us = usolution
        else: 
            
            #avalanche didn't occur, so find another and then go back
           
            print "avalanche didn't occur, trying to find next site to push"
            testModel.increaseW(mwpush/testModel.m**2)
            testModel.updateAllMatrix() # because the forces need to be updated
            # need to use dot here instead of multiply!
            # supposed to be 
            #ut = lambda t: np.dot(np.dot(np.dot(scipy.linalg.inv(testModel.evects),np.diag(np.exp(testModel.evals*t))), testModel.evects), testModel.us)-testModel.rhs*t
            ut = lambda t: solve(testModel.evects, np.diag(np.exp(testModel.evals*t))).dot(testModel.evects).dot(testModel.us) - testModel.rhs*t  
            # evaluate ut a little bit and see who runs into u_next first?
            # integrator not working very well yet... 
            t = 0.0
            increment_t = 1.0
            sols = testModel.us
            while (sols < u_upper).all():
                previous_t = t
                t += increment_t
                sols = ut(t)
                if np.sum(sols >= u_upper) > 1:
                    sols = testModel.us
                    increment_t = increment_t/2.
                    t = previous_t
            # now we've found which 
            minIndex = np.where(sols >= u_upper)[0][0]            
            #mwpush, minIndex = testModel.findFirstMinPush()
            #testModel.us[minIndex] = np.floor(testModel.us[minIndex])+1
            #u_next = np.fix(testModel.us)+1
            testModel.updatePartialMatrix([minIndex])
    
    print "now increase w"
    testModel.increaseW(mwpush)
    print "solution now:", testModel.us 

    print testModel.us
    print mwpush, minIndex
    while newpush <= mwpush:
        listmoved.append(newmin)
        print newmin, newpush
        testModel.us[newmin] = np.floor(testModel.us[newmin]) + 1
        newpush, newmin = testModel.findMinPush()
    print testModel.us   
    print listmoved
    testModel.increaseW(mwpush)
    testModel.updateAllMatrix()
    #testModel.updatePartialMatrix(listmoved)
    usolution = testModel.solveForUs()
    
    if (usolution < np.floor(testModel.us)).any() or (usolution >np.ceil(testModel.us)).any():
        print "solution not within range"
        # perhaps I should throw an exception here
    else:
        testModel.us = usolution
"""
