from scipy import integrate
from numpy import exp, inf, pi
class tests:
    def __init__(self,ep,p):
        self.ep = ep
        self.n = len(ep)
        self.p = p

    def CHF(self,a,b,z,*mode): #confluent Hypergeometric Functionm
        main = self.gamma(b)/(self.gamma(a)*self.gamma(b-a)) #the part out side the integral
        inte = integrate.quad(lambda t: exp(z*t)*(t**(a-1))*((1-t)**(b-a-1)),0,1)[0] #integrates the function between 1 and 0 giving numerical answer
        out = main*inte
        if mode[0] == "test":
            return out,main,inte;
        return out;

    def gamma(z): #donoted upper case gamma in spec
        Gam = integrate.quad(lambda t: (t**(z-1))*(exp(-t)),0,inf)[0]
        return Gam;

    def igame(self,mode,a,x): #incompete gamma function
        mode = mode.upper() #converts mode to upper case to reduce stupid bugs
        gam = self.Gamma(a)
        if mode == "Q":
            upper = x
            lower = 0
        elif mode == "P":
            upper = inf
            lower = x
        inte = integrate.quad(lambda t: exp(-t)*(t**(a-1)),lower,upper)

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
        P = igame("Q",N/2,chis/2)
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
		ep = 1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000
		Test = tests(ep,0.1)
		out = Test.BlockFrequency(10)
		if round(out,6) != 0.706438:
  			raise ValueError("Incorrect output BlockFrequency test 1 given{0}, needed{1}".format(out,0.706438))
  	
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
test_Frequency()
test_inte()
test_erfc()
