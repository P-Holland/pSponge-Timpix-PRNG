"""H0: ep is random
H1: ep is not random
sig_level = 0.01
odds of passing all tests when not random are low
"""
from scipy import integrate
from numpy import exp, inf, pi
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

    def gamma(self,z): #donoted upper case gamma in spec
        Gam = integrate.quad(lambda t: (t**(z-1))*(exp(-t)),0,inf)[0]
        return Gam;

    def igame(self,mode,a,x): #incompete gamma function
        mode = mode.upper() #converts mode to upper case to reduce stupid bugs
        gam = self.gamma(a)
        if mode == "P":
            upper = x
            lower = 0
        elif mode == "Q":
            upper = inf
            lower = x
        inte = integrate.quad(lambda t: exp(-t)*(t**(a-1)),lower,upper)
        return inte[0]/gam;

    def erfc(self,a):
        #print(others)
        out = 2/(pi**(1/2))
        in_ = integrate.quad(lambda u: exp(-(u**2)),a,inf)[0]
        whole = out*in_
        return whole;

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
        Pi  = []
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

    def LongestRunOfOnes(self,M):
        """M of 8,128 or 10000 are recommended, K of 3,5 or 6 are recommended (match to the similar M ie 8 and 3)"""
        #n = MN, N = n/M
        ep = self.ep
        n = self.n
        if n<750000:
            print("A longer string is needed for accurate results, or a low M value")
        N = n/M
        N = N-(N%1) #gets the floor of n/M
        ep = ep[:N*M] #shortens ep to the needed length making it easier to use (less proccessing required)
        a = []
        for i in range(M+1): #creates an empty list to store values
            a.append(0)
        for i in range(N):
            current = 0
            highest = 0
            in_use = ep[M*i:M*(i+1)]
            for b in in_use:
                if b == "1":
                    current +=1
                else:
                    if current>highest:
                        highest = current
                    current = 0
            a[highest]+=1
        if M == 8:
            v_map = {1:0,2:1,3:2,4:3}
            v = [  0, 0 , 0 , 0 ]
            low,high = 1,4
            Pi = [0.2148,0.3672,0.2305,0.1875]
            K = 3
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
        elif M = 10000:
            v_map = {10:0,11:1,12:2,13:3,14:4,15:5,16:6}
            v = [  0 ,  0 ,  0 , 0  , 0  ,  0 , 0  ]
            low,high = 10,16
            Pi = [0.0882,0.2092,0.2483,0.1933,0.1208,0.0675,0.0727]
            K = 6
        for i in range(0,K+1):
            if i<low:
                v[0]+=a[i] 
            elif i>high:
                v[len-1]+=a[i]
            else:
                v[v_map[i]]+=a[i] 
        chis = 0
        for i in range(K+1):
            chis+=(((v[i]-(N*Pi[i]))**2)/(N*Pi[i]))
        P = self.igame(K/2,chis/2)
        return P;


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

test_Frequency()
test_inte()
test_erfc()
test_blockFrequency()
test_runs()
