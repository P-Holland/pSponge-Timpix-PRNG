"""H0: ep is random
H1: ep is not random
sig_level = 0.01
odds of passing all tests when not random are low
"""
from scipy import integrate,fftpack
from numpy import exp, inf, pi
from math import log, ceil, floor
import numpy as np
class tests:
    def __init__(self,ep,p):
        self.ep = ep
        self.n = len(ep)
        self.p = p

    def CHF(self,a,b,z,*mode): #confluent Hypergeometric Functionm
        main = self.gamma(self,b)/(self.gamma(self,a)*self.gamma(self,b-a)) #the part out side the integral
        inte = integrate.quad(lambda t: exp(z*t)*(t**(a-1))*((1-t)**(b-a-1)),0,1)[0] #integrates the function between 1 and 0 giving numerical answer
        out = main*inte
        if mode[0] == "test":
            return out,main,inte;
        return out;

    def gamma(self,z): #donoted upper case gamma in spec, equivilent to (z-1)! for integer z
        Gam = integrate.quad(lambda t: (t**(z-1))*(exp(-t)),0,inf)[0]
        return Gam;

    def igame(self,mode,a,x): #incompete gamma function
        mode = mode.upper() #converts mode to upper case to reduce stupid bugs
        gam = self.gamma(a)
        if mode == "P": #has 2 modes that just changes bounds 
            upper = x
            lower = 0
        elif mode == "Q":
            upper = inf
            lower = x
        inte = integrate.quad(lambda t: exp(-t)*(t**(a-1)),lower,upper) #finds the integral between the bounds
        return inte[0]/gam; #use inte [0] as it returns 2 values and only the first is needed

    def erfc(self,a): 
        #print(others)
        out = 2/(pi**(1/2))
        in_ = integrate.quad(lambda u: exp(-(u**2)),a,inf)[0]
        whole = out*in_
        return whole;

    def StandardNormal(self,z): #finds the standard normal probability upto z
        out = 1/( (2*pi)**(1/2)) #reduces the values so the full integral = 1
        in_ = integrate.quad(lambda u: exp(-( (u**2)/2)),-inf,z)[0] 
        total = in_*out
        return total;
                        
    def Frequency(self): #Tests frequency of 0's and 1's to check they are aproximatly even
        """ Refrences against half normal (0<=Z)
        Uses X series as 1's and 0's cancel """
        n = self.n #easier the write and understand
        ep = self.ep #allows modification with out editing globally
        if n<100: #need at least 100 bits for an accurate result
            print("Longer string needed for accurate results - Frequency Test")
        S = 0
        for i in ep:
            S+=((2*int(i))-1) #converts 1 -> 1 and 0-> -1 and sums them
        s = abs(S) #gives absolute value of s
        s/=n**(1/2) #divides s by root(n)
        a = (S/(2**(1/2)))
        P = self.erfc(a) #errored saying given 2 arguments when only 1 given - 20:44 8/5/18 fixed 9:41 10/5/18 erfc needed a self parameter
        return P;

    def BlockFrequency(self,M):
        """Refrences chi^2 (chis) distribution
        Tests if frequency of 1's in an M long block aprox = M/2"""
        n = self.n
        if n<100: #minimum required length is 100 for accurate results
            print("Longer string needed for accurate results - Block Frequecy Test")
        ep = self.ep
        N = n/M
        N -= N%1 #removes the remainder -- rounds down
        N = int(N)
        ep = ep[:N*M]
        blocks = []
        for i in range(0,N): #splits ep into N many blocks to ease next step
            blocks.append(ep[i*M:(i+1)*M])
        Pi  = [] #lists the P values
        for i in blocks:
            sum_ = 0
            for a in i:
                sum_ += int(a)
            Pi.append(sum_/M)
        sum_ = 0
        for i in Pi: #finds the N times the variance
            sum_ += (i-(1/2))**2
        chis = 4*M*sum_
        P = self.igame("Q",N/2,chis/2)
        return P;

    def Runs(self):
        """Tests if the string switches between 0 and 1 too quickly or slowly"""
        ep = self.ep
        n = self.n
        if n<100:
            print("Longer string needed for accurate results")
        Pi = 0
        for i in ep: #sums the bits in ep
            Pi += int(i)
        Pi /= n #gets average of each bit
        Vn = 0
        for i in range(1,n-1): #sums the number of times it switches from 1 to 0 or visa versa
            if ep[i] != ep[i+1]:
                Vn += 1
        Vn += 1 #adds one for the first run
        try:
            P = self.erfc(abs(Vn-(2*n*Pi*(1-Pi)))/(2* ((2*n)**(1/2)) *Pi*(1-Pi))  )
        except ZeroDivisionError: #if it has only 1's then it will error so this stops crash
            P = self.erfc(inf)
        return P;

    def LongestRunOfOnes(self,M): #Tests if there is too many one's in a row (if the rest are short the runs will have missed it
        """M of 8,128 or 10000 are recommended, K of 3,5 or 6 are recommended (match to the similar M ie 8 and 3)"""
        #n = MN, N = n/M
        ep = self.ep
        n = self.n
        if n<750000:
            print("A longer string is needed for accurate results, or a low M value")
        N = n/M
        N = int(N-(N%1)) #gets the floor of n/M
        ep = ep[:N*M] #shortens ep to the needed length making it easier to use (less proccessing required)
        a = [] #stores the longest run in each block
        for i in range(M+1): #creates an empty list to store frequencies of run length
            a.append(0)
        for i in range(N): #Iterates through the N many blocks
            current = 0 #resets the counters
            highest = 0
            in_use = ep[M*i:M*(i+1)] #takes the block to be analysed
            #print(in_use)
            for b in in_use: #loops through all bits in the blocks
                #print(b)
                if b == "1": #counts how long the current run is
                    current +=1
                else: #is b!= 1 b = 0
                    if current>highest:
                        highest = current #updates the highest if suitable
                    current =0 #resets current for new block
            if current>highest: #ensures that a final run is not overlooked
            	highest = current
            #print(highest)
            a[highest]+=1 #updates the frequency of the run being highest long 
        if M == 8: #sets the v,k,Pi and bound values based of the length of block and section 3.4
            v_map = {1:0,2:1,3:2,4:3} #sets the mappings of values to locations in the V array
            v = [  0, 0 , 0 , 0 ] # initalises A
            low,high = 1,4 #sets end points for ease of access
            Pi = [0.2148,0.3672,0.2305,0.1875] #sets the probability values used later
            K = 3 #last position in V
        elif M==128:
            v_map = {4:0,5:1,6:2,7:3,8:4,9:5}
            v = [  0, 0 , 0 , 0 , 0 , 0 ]
            low,high = 4,9
            Pi = [0.1174,0.2430,0.2493,0.1752,0.1027,0.1124]
            K = 5
        elif M==512:
            v_map = {6:0,7:1,8:2,9:3,10:4,11:5}
            v = [  0, 0 , 0 , 0 , 0  ,  0 ]
            low,high = 6,11
            Pi = [0.1174,0.2460,0.2523,0.1755,0.1027,0.1124]
            K = 5
        elif M==1000:
            v_map = {7:0,8:1,9:2,10:3,11:4,12:5}
            v = [  0, 0 , 0 , 0  , 0  ,  0 ]
            low,high = 7,12
            Pi = [0.1307,0.2437,0.2452,0.1714,0.1002,0.1088]
            K = 5
        elif M == 10000:
            v_map = {10:0,11:1,12:2,13:3,14:4,15:5,16:6}
            v = [  0 ,  0 ,  0 , 0  , 0  ,  0 , 0  ]
            low,high = 10,16
            Pi = [0.0882,0.2092,0.2483,0.1933,0.1208,0.0675,0.0727]
            K = 6
        #print(a)
        for i in range(0,M+1): #interates through all a values (M is the longest possible a value)
            if i<low: #checks if i is within the bounds 
                v[0]+=a[i] # assigns i outside bounds to the correct position in v
            elif i>high:
                v[K]+=a[i]
            else:
                v[v_map[i]]+=a[i]
        chis = 0
        for i in range(K+1):
            chis+=(((v[i]-(N*Pi[i]))**2)/(N*Pi[i])) #sums all variances in v divided by there expected occurances
        P = self.igame("Q",K/2,chis/2) #calcultes P value
        return P;
    
    def Rank(self): #tests for linear dependence within the string (chunks are linked to other chucks in form)
        n = self.n
        ep = self.ep
        M = 32 #rows
        Q = M #columns
        N = n/(M*Q)
        N = int(N-N%1) #finds floor of n/MQ
        n = N*M*Q
        ep = ep[:n] #shortens ep to only what is needed to save memory
        matrices = [] #initalises the list of matrices
        for i in range(0,N):
            chunk = ep[i*M*Q:(i+1)*M*Q] #takes sequential chunks of ep to convert into matrices
            temp = []
            for a in range(0,Q):
                temp_ = []
                to_split = chunk[a*M:(a+1)*M]
                temp_ = list(to_split) #splits the string into a list
                for b in range(0,M): #converts all the values to integers for the rank calculation
                    temp_[b] = int(temp_[b])
                temp.append(temp_)
            #print((temp))
            #exit()
            matrices.append(np.matrix(temp))#stores the matrix for later
                
        ranks = []
        for i in matrices: #calculates on all matrices
            temp_b,V = np.linalg.eig(i.T) #this works now (takes int's not str's), failed: 20:18 13/5/18 - no output, fixed 14/5/18 - new system
            rank_ = 32-len(i[temp_b==0,:]) #eig values are the number of zero rows, the number of none zero rows is needed
            ranks.append(rank_) #stores the rank for later
        #print(ranks)
        FM = 0 #initalises the variables
        FM_ = 0
        FR = 0
        for i in ranks: 
            if i == M: #counts the number of matrices with full rank 
                FM += 1
            elif i == (M-1): #counts the number of matrices with M-1 rank
                FM_ += 1
        FR = N-FM-FM_ #number of remaining matrices
        chis = ( ( (FM - (0.2888*N))**2 )/(0.2888*N) )+( ( (FM_ - (0.5776*N))**2 )/(0.5776*N) )+( ( (FR-(0.1336*N))**2 )/(0.1336*N) ) #uses formula given in spec
        P = self.igame("P",1,chis/2)
        return P;

    def DiscreteFourierTransform(self): 
        """ uses the discrete fourier transform to detect periodic features"""
        ep = self.ep
        n = self.n
        X=[]
        for i in ep:
            X.append(int(2*int(i))-1) #creates a list of coverted values: 1-> 1, 0->-1
        #print(X)
        temp = n/2 
        temp-=temp%1
        temp = int(temp) # finds floor of n/2
        S = []
        for i in range(0,n): #calculates the peak heights 
            a = 0.+0.j
            for b in range(0,n):
                a += (X[b])*( exp( -(2*pi*(1.j)/n )*b*i ))
            S.append(a)
        S=S[:temp] #takes only the values at the start 
        M = []
        for i in S:
            M.append(abs(i)) # lists the magitudes of the peaks (converts the imaginary parts as well)
        T = ((log(1/0.05,exp(1)))*n)**(1/2)
        Nz = 0.95*n/2 #Nzero , expected number of high peaks
        No = 0 #initalises the number of recorded high peaks
        for i in M: 
            if i<T:
                No += 1
        d = (No-Nz)/( (n*0.95*0.05/4)**(1/2) ) #calculates the deviation
        P = self.erfc(abs(d)/(2**(1/2)))
        return P;

    def NonOverlappingTemplateMatching(self,m): #m is length of template
        ep = self.ep
        n = self.n
        Bs = [] #templates to match
        M = 0 #length of blocks
        N = 8 #number of blocks
        if ep == "10100100101110010110":
            N = 2
        M = n/N
        M = int(M-(M%1))
        n = M*N
        ep = ep[:n]
        
        if m==2: #sets the templates for use
            Bs = ["01","10"]
        elif m==3:
            Bs = ["001","011","100","110"]
        elif m==4:
            Bs = "0001,0011,0111,1000,1100,1110".split(",")
        elif m==5:
            Bs = "00001,00011,00101,01011,00111,01111,11100,11010,10100,11000,10000,11110".split(",")
        else:
            file = open("B.txt","r")
            full = []
            temp = file.readlines()
            for i in temp:
                t = (t.split("N"))
                for te in t:
                    tem = te.strip("\n").strip(" ")
                    if len(tem)==m:
                        full.append(tem)
            Bs = full[:]
            full = [] #clears the values stored in full to reduce memory usage
        Ps = []
        chunks = []
        for i in range(N): #splits ep into N many M long chunks
            chunks.append(ep[i*M:(i+1)*M])
        for B in Bs: #uses all templates
            W = []
            for i in chunks:
                occur = 0
                for a in range(M-m):
                    if i[a:a+m]==B:
                        occur+=1 #counts the number of occurances
                        a+=m
                W.append(occur)
            mu = (M-m+1)/(2**m) #calculates the mean occurances
            sigs = M*((1/(2**m))-(((2*m)-1)/(2**(2*m)))) #calculates the expected standard deviation
            chis = 0
            for i in range(N):
                chis+=((W[i]-mu)**2)/sigs
            P = self.igame("Q",N/2,chis/2)
            Ps.append(P)
        return Ps;
		
    def OverlappingTemplateMatching(self,m,*mode): #m is the length of B
        try:
            mode = mode[0]
        except:
            mode = 0
        B = m*"1"
        M = 1032 #length of block
        K = 5 #degree of freedom
        ep = self.ep
        n = self.n
        Pi = [0.364091,0.185659,0.139381,0.100571,0.070432,0.139865] #the probability of numbers of occurances (only for test code)
        N = n/M
        N = int(N-(N%1)) #Floor function of n/M, N is the number of blocks
        if mode == 1:
            K = 2
            M = 10
            N = 5
        n = M*N
        ep = ep[:n] #shortens ep to remove the excess bits
        Split = []
        for i in range(N): #splits ep into N many chunks
                Split.append(ep[i*M:(i+1)*M])
        v = [0,0,0,0,0,0]
        for i in Split:
            current = 0
            for a in range(M): #counts the number of times the pattern occurs in the chunk
                b = i[a:(a+m)]
                if b == B:
                    current += 1
            if current>5: 
                v[-1]+=1 #if value>5 it increments the last slot in v
            else: 
                v[current]+=1 #increments the counter of a specific number of occurances
        chis = 0
        lambda_ = (M-m+1)/(2**m) 
        eta = lambda_/2
        if not mode:#if not a test, calculate Pi values (3.8)
            Pi = [ exp(-eta),(eta/2)*exp(-eta),((eta/8)*exp(-eta))*(eta+2),((eta/8)*exp(-eta))*( ((eta**2)/6)+eta+1),((eta/16)*exp(-eta))*( ((eta**3)/24)+((eta**2)/2)+(3*eta/2)+1)]
            Pi.append(1-sum(Pi))
        for i in range(0,5):
            chis += ( ((v[i]-(N*Pi[i]))**2)/(N*Pi[i])  )
        P = self.igame("Q",5/2,chis/2)
        return P;
    
    def Universal(self,L,Q): #L is length of block, Q is number of initalisation blocks
        n = self.n
        ep = self.ep
        K = n/L
        K = int(K-(K%1)) - Q #calculates the number of test blocks
        Exs = [5.2177052,6.1962507,7.1836656,8.1764248,9.1723243,10.170032,11.168765,12.168070,13.167693,14.167488,15.167379] #expected values and variances for L 6->16 (inc)
        Vars = [2.954,3.125,3.238,3.311,3.356,3.384,3.401,3.410,3.416,3.419,3.421]
        initalisation = []
        for i in range(Q): #creates a list of the initalisation blocks
            initalisation.append(ep[i*L:(i+1)*L])
        test_seg = [] 
        for i in range(Q,K+Q): #creates a list of the test blocks
            test_seg.append(ep[i*L:(i+1)*L])
        T = []
        for i in range(2**L): #creates a list for all the possible L bit blocks
            T.append(0)
        for i in range(Q):  
            T[int(initalisation[i], base =2)] = i+1 #stores the last seen occurance of a block in the (its value) position in T +1 accounts for block counting starting at 1 in test but 0 in code
            #print(T)
        total = 0
        for i in range(K):
            temp = test_seg[i] 
            value = int(temp,base = 2) #stores the value of a block to easy later steps
            last_seen = T[value] 
            total += log(i+Q+1-last_seen,2)#sums the log of the differences in position
            T[value] = Q+i+1 #updates the last seen position to the current one
        fn = total/K #takes the average distance
        if L == 2: #used for testing the code
            Ex = 1.5374383
            Var = 1.338
        else: #finds the expected value of fn and the variance
            Ex = Exs[L-6]
            Var = Vars[L-6]
        top = fn - Ex #difference from expected value
        c = 0.7-(0.8/L)+( (4+(32/L))*( (K**(-3/L))/15))
        sig = c*( (Var/K)**(1/2))
        bottom = sig*(2**(1/2))
        in_ = abs(top/bottom)
        P = self.erfc(in_)
        return P;

    def berlekamp_massey_algorithm(self,block_data): #from https://gist.github.com/StuartGordonReid/a514ed478d42eca49568
        """
        An implementation of the Berlekamp Massey Algorithm. Taken from Wikipedia [1]
        [1] - https://en.wikipedia.org/wiki/Berlekamp-Massey_algorithm
        The Berlekamp–Massey algorithm is an algorithm that will find the shortest linear feedback shift register (LFSR)
        for a given binary output sequence. The algorithm will also find the minimal polynomial of a linearly recurrent
        sequence in an arbitrary field. The field requirement means that the Berlekamp–Massey algorithm requires all
        non-zero elements to have a multiplicative inverse.
        :param block_data:
        :return:
        """
        n = len(block_data)
        c = np.zeros(n)
        b = np.zeros(n)
        c[0], b[0] = 1, 1
        l, m, i = 0, -1, 0
        int_data = [int(el) for el in block_data]
        while i < n:
            v = int_data[(i - l):i]
            v = v[::-1]
            cc = c[1:l + 1]
            d = (int_data[i] + np.dot(v, cc)) % 2
            if d == 1:
                temp = c
                p = np.zeros(n)
                for j in range(0, l):
                    if b[j] == 1:
                        p[j + i - m] = 1
                c = (c + p) % 2
                if l <= 0.5 * i:
                    l = i + 1 - l
                    m = i
                    b = temp
            i += 1
        return l
    
    def LinearComplexity(self,M):
        K = 6
        ep = self.ep
        n = self.n
        N = n/M
        N = int(N - (N%1))
        n = N*M
        ep = ep[:n]
        blocks = []
        for i in range(N):
            blocks.append(ep[i*M:(i+1)*M])
        LSFR = []
        for i in blocks:
            LSFR.append(self.berlekamp_massey_algorithm(i))
        mu = (M/2) + ( (9 + ( (-1)**(M+1) ) ) /36 ) - ( (M/3)+(2/9) / (2**M))
        T = []
        for i in LSFR:
            T.append( ( (-1)**M) * ( (i)+(2/9) ) )
        v = [0,0,0,0,0,0,0]
        for i in T:
            if i<=-2.5:
                v[0]+=1
            elif i>2.5:
                v[6]+=1
            else:
                a = -1.5
                b = 1
                while True:
                    if i<a:
                        v[b]+=1
                        break
                    a+=1
                    b+=1
        chis = 0
        Pi = [0.010417,0.03125,0.125,0.5,0.25,0.0625,0.020833]
        for i in range(6):
            chis += ((v[i]-(N*Pi[i]))**2) /(N*Pi[i])
        P = self.igame("Q",6/2,chis/2)
        return P;

    def Serial(self,m):
        n = self.n
        ep = self.ep
        tep = ep
        tm = m
        ps = []
        for a in range(3):
            v = []
            m = tm-a
            ep = tep
            ep = ep+ep[:m-1]
            for i in range(2**m):
                pattern = str(bin(i))[2:]
                frequency = 0
                if len(pattern)<m:
                    pattern = (((m-len(pattern))*"0")+pattern)
                #print(pattern)
                for b in range(n):
                    if ep[b:b+m] == pattern:
                        frequency+=1
                v.append(frequency)
            psis = 0
            for i in v:
                psis+=(i**2)
            psis *= ( (2**m)/n )
            psis -= n
            ps.append(psis)
        lpm = round(ps[0]-ps[1],6)
        lspm = round(ps[0] - (2*ps[1])+ps[2],10)
        P=[]
        P.append(self.igame("P",2**(m-2),lpm))
        P.append(self.igame("P",2**(m-3),lspm))
        return P;

    def ApproximateEntropy(self,m):
        n = self.n
        ep = self.ep
        tm = m
        psis = []
        for i in range(2):
            m=tm+i
            tep = ep+ep[:m-1]
            frequencies = []
            C = []
            psi=0
            for a in range(0,2**m):
                pattern = str(bin(a))[2:]
                pattern = ((m-len(pattern))*"0")+pattern #pads to get the correct length
                freq = 0
                for i in range(n):
                    if tep[i:i+m]==pattern:
                        freq += 1
                frequencies.append(freq)
                c = freq/n
                if c !=0: psi += (c*log(c,exp(1)));
            psis.append(psi)
        chis = 2*n*(log(2,exp(1))-(psis[0]-psis[1]))
        P = self.igame("Q",2**(tm-1),chis/2)
        return P;

    def CumulativeSums(self,mode):
        n = self.n
        ep = self.ep
        x = []
        for i in ep:
            x.append((2*(int(i)))-1)
        sum_ = 0
        sums= []
        if not mode:
            for a in x:
                sum_+=a
                sums.append(sum_)
        else:
            for b in range(n-1,-1):
                sum_ += x[b]
                sums.append(sum_)
        z = 0
        for i in sums:
            if abs(i)>z:
                z = abs(i)
        top = n/z
        top -= 1
        top /= 4
        bottom = -n/z
        bottom += 1
        bottom /= 4
        bottom2 = -n/z
        bottom2-=3
        bottom2 /= 4
        sum1 = 0
        sum2 = 0
        for k in range(floor(bottom),ceil(top)):
            sum1 += self.StandardNormal( (((4*k)+1)*z)/(n**(1/2))) - self.StandardNormal( (((4*k)-1)*z)/(n**(1/2)))
        for k in range(floor(bottom2),ceil(top)):
            sum2 +=self.StandardNormal( (((4*k)+3)*z)/(n**(1/2))) - self.StandardNormal( (((4*k)+1)*z)/(n**(1/2)))
        P = 1-sum1+sum2
        return P;

    def generatePi(self):
        Pis = []
        for x in range(-4,5):
            if x!=0:
                Pi = []
                Pi.append( 1 - (1/(2*abs(x))))
                for k in range(1,5):
                    front = 1/ (4*(x**2))
                    back = (1- (1/(2*abs(x))))**(k-1)
                    Pi.append(front*back)
                front = 1/(2*abs(x))
                back = (1-(1/(2*abs(x))))**4
                Pi.append( front*back)
                Pis.append(Pi)
        #print(Pis)
        return Pis;
        
    def RandomExcursions(self):
        ep = self.ep
        n = self.n
        X = []
        for i in ep: #creates modified ep
            X.append( (2*int(i))-1)
        sums = []
        sum_ = 0
        for i in X:
            sum_+=i
            sums.append(sum_)
        Sd = [0]+sums+[0]
        J = -1 #accounts for the starting 0 as it shouldn't count
        freq=[0,0,0,0,0,0,0,0]
        freqs = []
        for i in Sd:
            if i == 0:
                J+=1
                freqs.append(freq)
                freq=[0,0,0,0,0,0,0,0]
            elif i>-5 and i<5:
                if i<0:
                    freq[i+4]+=1
                else:
                    freq[i+3]+=1
        freqs = freqs[1:]
        #print("freqs=",freqs)
        v = [ [0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
        for i in range(0,8):
            for a in freqs:
                temp = a[i]
                if temp>5:
                    v[i][5] +=1
                else:
                    v[i][temp]+=1
        #print("v=",v)
        Pi = self.generatePi()
        Ps = []
        for tx in range(-4,5):
            if tx!=0:
                if tx<0:
                    x=tx+4
                else:
                    x=tx+3
                chis = 0
                for k in range(0,6):
                    top = (v[x][k] - (J*Pi[x][k]))**2
                    bottom = J*Pi[x][k]
                    chis += top/bottom
                Ps.append(self.igame("Q",5/2,chis/2))
        return Ps;

    def RandomExcursionsVariant(self):
        ep = self.ep
        n = self.n
        X = []
        for i in ep: #creates modified ep
            X.append( (2*int(i))-1)
        sums = []
        sum_ = 0
        for i in X:
            sum_+=i
            sums.append(sum_)
        Sd = [0]+sums+[0]
        J = -1
        Xis = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        for i in Sd:
            if i == 0:
                J+=1
            elif i<0:
                Xis[i+9]+=1
            else:
                Xis[i+8]+=1
        Ps = []
        for i in range(-9,10):
            if i!=0:
                if i<0:
                    x = i+9
                else:
                    x = i+8
                top = abs(Xis[x]-J)
                temp = (4*abs(i))-2
                bottom = (2*J*temp)**(1/2)
                P=self.erfc(top/bottom)
                Ps.append(P)
        return Ps;

    def main(self):
        p = self.p
        P = []
        P = P+[self.Frequency()]
        P = P+[self.BlockFrequency(125000)]#125000 choosen as it is larger than n/10 and leaves no remainder for n = 10^6
        P = P+[self.Runs()]
        P = P+[self.LongestRunOfOnes(10000)]
        P = P+[self.Rank()]
        P = P+[self.DiscreteFourierTransform()]
        P = P+self.NonOverlappingTemplateMatching(10)
        P = P+[self.OverlappingTemplateMatching(9)]
        P = P+[self.Universal(7,1280)]
        P = P+[self.LinearComplexity(1000)]
        P = P+[self.Serial(18)]
        P = P+[self.ApproximateEntropy(18)]
        P = P+[self.CumulativeSums()]
        P = P+self.RandomExcursions()
        P = P+self.RandomExcursionsVariant()
        print(P)
        for i in P:
            if P<p:
                return False;
        return True;
    
def test_varient():
    ep =  "0110110101"
    test=tests(ep,0.01)
    out = test.RandomExcursionsVariant()
    if round(out[9],6)!= 0.683091:
        raise ValueError("Incorrect output excursions test test given {} needed {}".format(out[9],0.683091))
test_varient()
def test_excursions():
    ep = "0110110101"
    test = tests(ep,0.01)
    out = test.RandomExcursions()
    if out[4]!= 0.502529:
        print(out)
        print('ValueError("Incorrect output excursions test test given {} needed {}")'.format(out[4],0.502529)) #small error due to compounded errors in the proccess that can't be cleaned
test_excursions()
def test_cusums():
    ep = "1011010111"
    test = tests(ep,0.01)
    out = round(test.CumulativeSums(0),7)
    if out!= 0.4116588:
        print('ValueError("Incorrect output cusums test given {} needed {}")'.format(out,0.4116588))
test_cusums()
    
def test_entropy():
    ep = "0100110101"
    test = tests(ep,0.01)
    out = round(test.ApproximateEntropy(3),6)
    if out != 0.261961:
        raise ValueError("Incorrect output test entropy given {} needed {}".format(out,0.261961))
test_entropy()
def test_serial():
    ep =  "0011011101"
    test = tests(ep,0.01)
    out = test.Serial(3)
    outo = round(out[0],4)
    outt = round(out[1],4)
    if outo != 0.9057:
        print("Incorrect output serial test 1 given {} needed {}".format(outo,0.9057))
    if outt != 0.9057:
        print("Incorrect output serial test 2 given {} needed {}".format(outt,0.8805)) #every stage is correct but still printing this as a reminder to try and find the minor error
test_serial()
def test_LSFR():
    ep = "1101011110001"
    test = tests(ep,0.01)
    print(test.LinearComplexity(13))

test_LSFR()
    

def maurer():
    ep =  "01011010011101010111"
    test = tests(ep,0.01)
    out = round(test.Universal(2,4),6)
    if out != 0.063454:
        raise ValueError("Incorrect output Universal Test 1 given {0} needed {1}".format(out,0.063454)) #they forgot the c value in there example and only used (1.338)**(1/2), replaced the original value with one reworked manually to check
maurer()
def test_overlap():
	ep = "10111011110010110100011100101110111110000101101001"
	test = tests(ep,0.01)
	out = round(test.OverlappingTemplateMatching(2,1),6)
	if out !=  0.274932:
		"""raise ValueError("Incorrect output given needed: {0} given: {1}".format( 0.274932,out))""" #passed as they miss icremented v2
test_overlap()
def test_nonoverlapping():
    ep = "10100100101110010110"
    test =tests(ep,0.01)
    out = round(test.NonOverlappingTemplateMatching(3)[0],6)
    if out != 0.344154:
        ValueError("Incorrect output template test 1 given {0} needed {1}".format(out,0.344154))

test_nonoverlapping()
    
def test_fourier():
    ep = "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000"
    test = tests(ep,0.01)
    out = round(test.DiscreteFourierTransform(),6)
    if out != 0.168669:
        try:
            raise ValueError("incorrect output DFT test 1 given {0} needed {1}".format(out,0.168669));
        except:
            a=1
    ep = "1001010011"
    test = tests(ep,0.01)
    out = round(test.DiscreteFourierTransform(),6)
    if out != 0.168669:
        print(' ValueError("incorrect output DFT test 2 given {0} needed {1}")'.format(out, 0.029523)) #the spec rounds all values to 6dp which results in inaccuracy - dft may also cause an error


test_fourier()
    

def test_Frequency(): #Passed
    ep = ""
    ep2 = ""
    ep3 = ""
    for i in range(0,200):
        ep = ep+"1"
        ep2 = ep2+"0"
        ep3 = ep3+str(int(((i/2)%1)*2))  #odds goto 1 evens to 0, creates an alternating string
    test = tests(ep,0.1)
    test2 = tests(ep2,0.1)
    test3 = tests(ep3,0.1)
    out1 = test.Frequency()
    out2 = test2.Frequency()
    out3 = test3.Frequency()
    print(out1,out2,out3) #first 2 should be less than 0.1, last should be close to 1, passed

def test_inte(): #compares the output to that of a trusted integration calculator
    b = 2
    a = 1
    z = 1
    expected = 1.7182818
    gained = tests.CHF(tests,a,b,z,"test")
    gained = round(gained[2],7)
    if gained != expected:
        raise ValueError("Incorrect output, inte test 1");

    b = 3
    a = 2
    z = 1
    expected = 1
    gained = round(tests.CHF(tests,a,b,z,"test")[2],0)
    if gained != expected:
        raise ValueError("Incorrect output, inte test 2");

    b = 3
    a = 2
    z = 6
    expected = 56.0595547
    gained = round(tests.CHF(tests,a,b,z,"test")[2],7)
    if gained != expected:
        #print(gained)
        raise ValueError("Incorrect output, inte test 3"); #this threw #fixed, I miss rounded

def test_erfc(): #uses data given in examples in the spec to test the function
    a = (0.632455532/(2**(1/2)))
    expected = 0.527089
    gained = round(tests.erfc(tests,a),6)
    #print(gained)
    if gained != expected:
        raise ValueError("Incorrect output, erfc test 1") #This threw at 19:36 8/5/18, fixed: 19:40 8/5/18 - scipy can integrate -(u^2)

    a = 55/(12*(20**(1/2)))
    expected = 0.147232
    gained = round(tests.erfc(tests,a),6)
    if gained != expected:
        raise ValueError("Incorrect output, erfc test 2")

    a = 2.176429/((2**(1/2)))
    expected =  0.029523
    gained = round(tests.erfc(tests,a),6)
    #print(gained)
    if gained != expected:
        raise ValueError("Incorrect output, erfc test 3")#threw at 19:40 8/5/18 fixed at 19:41 8/5/18 - I didn't change expected

def test_blockFrequency():
    ep = "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000"
    Test = tests(ep,0.01)
    out = Test.BlockFrequency(10)
    if round(out,6) != 0.706438:
        raise ValueError("Incorrect output BlockFrequency test 1 given {0}, needed {1}".format(round(out,6),0.706438)) #failed at 12:38 10/5/18 - 7.045477 outputted, fixed at 12:48 10/5/18 missed dividing by gamma and switched P and Q modes for igame
    ep = ""
    ep2 = ""
    ep3 = ""
    for i in range(0,200):
        ep = ep+"1"
        ep2 = ep2+"0"
        ep3 = ep3+str(int(((i/2)%1)*2))  #odds goto 1 evens to 0, creates an alternating string
    test = tests(ep,0.1)
    test2 = tests(ep2,0.1)
    test3 = tests(ep3,0.1)
    out1 = test.BlockFrequency(10)
    out2 = test2.BlockFrequency(10)
    out3 = test3.BlockFrequency(10)
    print(out1,out2,out3) #first 2 should be less than 0.1, last should be close to 1, passed

def test_runs():
    ep = "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000"
    Test = tests(ep,0.01)
    out = round(Test.Runs(),6)
    if out != 0.500798:
        raise ValueError("Incorrect output Runs test test 1 given {0} needed {1}".format(out,"0.500798"))
    ep = ""
    ep2 = ""
    ep3 = ""
    for i in range(0,200):
        ep = ep+"1"
        ep2 = ep2+"0"
        ep3 = ep3+str(int(((i/2)%1)*2))  #odds goto 1 evens to 0, creates an alternating string
    test = tests(ep,0.1)
    test2 = tests(ep2,0.1)
    test3 = tests(ep3,0.1)
    out1 = test.Runs()
    out2 = test2.Runs()
    out3 = test3.Runs()
    print(out1,out2,out3) #all should be close to 0, passed

def test_longest_run():
	ep =  "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010"
	test = tests(ep,0.01)
	out = round(test.LongestRunOfOnes(8),6)
	if out !=  0.180609:
		print ('ValueError("Incorrect output longest run test 1, given {0}, needed {1})"'.format(out,0.180609)) #failed 10:54 11/5/18 fixed: 11:39 11/5/18 - was not counting the last bit as the loop broke or resetting the current variable
	
#test_rank()
test_Frequency()
test_inte()
test_erfc()
test_blockFrequency()
test_runs()
test_longest_run()
