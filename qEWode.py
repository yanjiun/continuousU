import scipy.integrate
import numpy as np
import scipy.ndimage
import scipy.optimize
import scipy.sparse
import pylab
import pyfftw
from scipy.sparse.linalg import spsolve
from scipy.linalg import eig, solve
from scipy.integrate import odeint

# 
# 
# later on: find zeros of velocity with matrix method?
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

    def notify(self, exclude=None):
        for observer in self._observers:
            if observer not in exclude: # provides option to exclude some observers from being notified
                observer.update(self)


class qEWContinuous(Subject):
    def __init__(self,m,N,gamma,sigma,seed):
        Subject.__init__(self)
        np.random.seed(seed)
        self.randmat = sigma*np.random.randn(10*N,N)
        self.us = np.zeros(N)
        self.m = m
        self.w = 0.0
        self.length = N
        self.VelocityArray = np.zeros(N)
        self.lhs = np.zeros((N,N)) 
        self.gamma = gamma
        self.current_cs = 0
        self.current_as = 0    
        self.finish_run = False
        self.u_traj = 0

    def buildMatrix(self):
        """
        this will only be called once per simulation?
        do i build this directly as a sparse matrix or as a dense matrix then convert
        make a self.lhs and self.rhs? 
        """
        gammaRow = np.zeros(self.length) 
        gammaRow[1] = self.gamma
        gammaRow[-1] = self.gamma
        a = gammaRow
        for i in range(0,self.length-1):
            gammaRow = np.roll(gammaRow,1)
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
        # this is assuming that the current configuration is stable?
        mwpushes = self.m**2*(u_next-self.us) + 2*self.gamma*(u_next-self.us) -self.current_cs*(u_next-self.us)
        #maskedpushes = np.ma.array(mwpushes,mask=(mwpushes>0))
        #print maskedpushes
        print "local forces:", mwpushes
        minIndex = mwpushes.argmin()
        mwmin = mwpushes[minIndex]
            
        return mwmin, minIndex

    def findW(self,u0=None):
        
        if u0 is None:
            # this solves the current
            u0 = self.solveForUs()
            print u0
            
        deltau = spsolve(self.lhs,-np.ones(self.length))
        print "from find W deltau= ", deltau        

        deltaw = (np.ceil(u0)-u0)/(self.m**2*deltau)
        print "from find W detlaw= ", deltaw        

        return min(deltaw), deltaw.argmin()

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
        return -self.VelocityArray.min()/self.m**2.

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


    def calculateVelocityArray(self, u_array, tdummy=0.0):

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
        u_uppers = np.fix(self.us)+1 # use np.ceil
        
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

class storeVelocities:
    """
    store velocities of avalanche?
    """
    
    def __init__(self, fileName = "Velocities.txt"):
        self.fileObject = open(fileName, "w")

    def update(self, subject, sampling_freq = 100):
        """
        calculates the average velocity of the front throughout the avalanche at various time points in the avalanche
        """
        if not self.finish_run:
            duration, width = np.shape(subject.u_traj)    
            t_list = np.arange(0,duration,sampling_freq)
            average_velocities = [np.average(testModel.calculateVelocityArray(subject.u_traj[t])) for t in t_list]
            self.store()
        else:
            self.fileObject.close()

    def store(self):
        """
        stores the velocity information
        """
        self.fileObject.write(str(average_velocities)+"\n")


class storeHeights:
    """
    store heights at the places where the front stops (so F = 0 configurations)
    and where configuration is unstable...
    """
    def __init__(self, fileName = "Heights.txt"):
        self.fileObject = open(fileName, "w")

    def update(self,subject):
        if self.finish_run:
            self.fileObject.close()        
        else:
            self.fileObject.write(str(subject.us)+"\n")

class storeCorrelations:
    """
    calculate height height correlations for each solution of the model
    """
    def __init__(self, subject, fileName="Correlations.pkl"):
        self.fileName = fileName
        self.fileObject = open(self.fileName, 'w')
        self.correlations = np.zeros(len(subject.us))
        self.old_heights = np.zeros(len(subject.us))
        self.front_counter = 0
        self.total_weight = 0

    def calculate(self,subject):
        # let's calculate the unweighted correlation for now...
        #weight = np.sum(subject.us-self.old_heights)
        #self.correlations += weight*np.fft.fft(subject.us)
        
        # below is pyfftw code to be used for speed up on ept
        #shape_tuple = subject.us.shape
        #a = pyfftw.n_byte_align_empty(shape_tuple,16,dtype=np.complex128)
        #b = pyfftw.n_byte_align_empty(shape_tuple,16,dtype=np.complex128)
        #c = pyfftw.n_byte_align_empty(shape_tuple,16,dtype=np.complex128)

        #fft = pyfftw.FFTW(a,b)
    
        #a[:] = subject.us
        #fft.execute()
    
        #ifft = pyfftw.FFTW(b*b.conjugate(),c, direction="FFTW_BACKWARD")

        #ifft.execute()

        us_fft = np.fft.fft(subject.us) 
        c = np.fft.ifft(us_fft*us_fft.conjugate())

        av_corr = (c-np.average(subject.us))/subject.length        

        self.correlations += av_corr
        # for some reason I used to do the negative of this... why???
        self.sum_squares += av_corr**2.       
 
        self.front_counter += 1
        #self.total_weight += weight

    def update(self,subject): 
        if not subject.finish_run:
            self.calculate(subject) # calculate the correlation each time...
        # think about how to store this part...        
        if subject.finish_run:
            self.correlations=self.correlations/self.front_counter
            self.sum_squares = self.sum_squares/self.front_counter
            error = np.sqrt((self.sum_squares-self.correlations**2.)/(self.front_counter-1.))
            pickle.dump(self.fileObject,[self.correlations, error])


def main(m,N,gamma,sigma,seed=1):
    
    # establish simulation class
    testModel = qEWContinuous(m,N,gamma,sigma,seed)

    # attach observers
    dir = 'data/'
    commonName = 'm='+str(m) + 'N='+str(N)+'g='+str(gamma)
    
    heightsFile = dir + 'qEWheights' + commonName +'.txt'
    heightsObserver = storeHeights(filename = heightsFile)
    velocitiesFile = dir + 'qEWvelocities' + commonName + '.txt'
    velocitiesObserver = storeVelocities(filename = velocitiesFile)
    correlationsFile = dir + 'qEWcorr'+commonName +'.bp'
    corrObserver = storeCorrelations(filename = correlationsFile)
    
    testModel.attach(heightsObserver)
    testModel.attach(velocitiesObserver)
    testModel.attach(corrObserver)

    # set initial external force equal to the most negative a_s
    # increase external force slowly and solve the differential equation (monitor velocity for when avalanche stops

    testModel.VelocityArray = testModel.calculateVelocityArray(testModel.us)
    firstW = testModel.findDeltaW() # this finds the most negative a_s

    #testModel.buildMatrix()
    # put external force equal to most negative velocity
    # not optimal -- try pushing less, or more intelligently...
    #firstW, argchanged = testModel.findW()
   
    #testModel.increaseW(firstW/2)

    # now start integrating forward until front stops (what time points should we use)
    # determine

    # the integrateFront function integrates outside of the class... 
    # conisder logic   
    testModel.u_traj, topReached = integrateFront(testModel, firstW/2, 10.0,1000) 
    #pylab.figure()
    #pylab.plot(u_traj.transpose())
    plotTrajectories(testModel.u_traj)
    #plotVelocities(testModel,u_traj)
    testModel.us = u_traj[-1]
    testModel.buildMatrix()
    #print "testModel.us from main: ", testModel.us
    #print testModel.doesAvalancheHappen()

    #print testModel.calculateVelocityArray(testModel.us)
    #print np.dot(testModel.lhs.todense(),testModel.us) - testModel.rhs   

    deltaw, argchanged = testModel.findW()

    testModel.increaseW(deltaw)
    testModel.rhs = testModel.rhs - testModel.m**2*deltaw
    # line above should be consolidated in increaseW? (does ODE solver depend on this?)
    testModel.us = testModel.solveForUs()
    testModel.us[argchanged] = np.ceil(testModel.us[argchanged])
    testModel.updatePartialMatrix([argchanged])

    #print testModel.doesAvalancheHappen()
    print topReached


    while not (testModel.us >= (testModel.length*10-1)).any() and topReached is False:
        # repeat this portion until end of the system...
        if testModel.doesAvalancheHappen():
            print "avalanche happened!"
            # if avalanche happens, notify observers
            testModel.notify() # observers can be skipped...
            # if avalanche happens, push front a little bit...
            u_traj_new, topReached = integrateFront(testModel,0.001/m**2,10.0,1000) 
            #plotTrajectories(u_traj_new)
            #plotVelocities(testModel,u_traj_new)
            if not topReached:
                testModel.us = u_traj_new[-1]
                testModel.u_traj = u_traj_new
                testModel.updateAllMatrix()
                testModel.us = testModel.solveForUs()
                testModel.notify(exclude = [heightsObserver, corrObserver])
            else:
                testModel.finish_run = True
                testModel.notify()
            #u_traj = np.concatenate((u_traj,u_traj_new))  
        else:
            print "front creeping"
            deltaw, argchanged = testModel.findW() 
            testModel.increaseW(deltaw)
            testModel.rhs = testModel.rhs - testModel.m**2*deltaw        
            testModel.us = testModel.solveForUs()
            testModel.us[argchanged] = np.ceil(testModel.us[argchanged])
            testModel.updatePartialMatrix([argchanged])
    
         
    
    #    while not testModel.doesAvalancheHappen():
    #        try:
    #            u_traj_new = integrateFront(testModel,0.01/m**2,1.0,100)
    #            #pylab.figure()
    #            #pylab.plot(u_traj_new.transpose())
    #            #plotTrajectories(u_traj_new)
    #            #plotVelocities(testModel,u_traj_new)
    #            testModel.us = u_traj_new[-1]
    #            testModel.updateAllMatrix()
    #            print testModel.us
    #            print testModel.solveForUs()
    #            print testModel.doesAvalancheHappen()
    #            u_traj = np.concatenate((u_traj,u_traj_new))
                # now increase W slowly and integrate slowly
    #        except:
    #            break   
     

    

    #return testModel, u_traj

def integrateFront(testModel, Wincrements, t_end,time_steps):
    testModel.increaseW(Wincrements) 
    tarray = np.linspace(0.0,t_end, time_steps)        
    u_traj = odeint(testModel.calculateVelocityArray,testModel.us,tarray)
    # NOTE: need more elegant way to stop integration once out of range...
    topReached = False
    while np.average(testModel.calculateVelocityArray(u_traj[-1])) > 10**(-5):
        print "front still going"
        try:
            u_traj_new = odeint(testModel.calculateVelocityArray,u_traj[-1],tarray) 
            u_traj=np.concatenate((u_traj,u_traj_new))
            if (u_traj[-1]>=(testModel.length*10-1)).any():
                print "from integrateFront last trajectory: ", u_traj[-1]
                topReached = True
                return u_traj, topReached
        except:
            print "exception occured"
            topReached = True
            return u_traj, topReached    

    return u_traj, topReached

# Note: should I make another class the interacts with both subject and observer to do plotting?

def plotTrajectories(u_traj):
    n = len(u_traj)
    ind = np.arange(0,int(n/10)*10,10)
    pylab.figure()
    pylab.plot(u_traj[ind,:].transpose())      

def plotVelocities(testModel, u_traj):
    """
    plot a random site's local velocity
    also average velocities
    """

    duration, width = np.shape(u_traj)    

    site_to_track = np.random.randint(0,10)

    t_list = np.arange(0,duration-1,10)

    velocities = [testModel.calculateVelocityArray(u_traj[t])[site_to_track] for t in t_list]    
    average_velocities = [np.average(testModel.calculateVelocityArray(u_traj[t])) for t in t_list]

    pylab.figure()
    pylab.plot(velocities, 'r', label='single site velocity')
    pylab.plot(average_velocities, 'b', label='average velocity of front')
    pylab.legend() 

    # add functionality to track velocities as a function of external force    
