from math import log
import gc
class Sponge:
    def __init__(self,b,nr=0):
        self.iter = 0
        self.b = b
        self.w = int(self.b/25)
        self.fl = int(log(self.w,2))
        if nr == 0:
            nr = 12+(2*self.fl)
        self.nr = nr
        
    def S_to_A(self,S):
        A  = self.form_Ad()
        w = self.w
        x = 0
        y = 0
        for x in range(5):
            for y in range(5):
                for z in range(w):
                    A[x][y][z]=S[(w*( (5*y)+x))+z]
        gc.collect()
        return A;
    
    def A_to_S(self,A):
        #print(len(A))
        d = []
        for a in range(0,5):
            d.append("")
        Lane = []
        for b in range(0,5):
            Lane.append(d[:])
        for i in range(0,5):
            for j in range(0,5):
                temp = ""
                for b in range(0,self.w):
                    temp = temp+str(A[i][j][b])
                Lane[i][j]=temp
        Plane = ["","","","",""]
        for j in range(0,5):
            Plane[j]=str(Lane[0][j])+str(Lane[1][j])+str(Lane[2][j])+str(Lane[3][j])+str(Lane[4][j])
        S = str(Plane[0])+str(Plane[1])+str(Plane[2])+str(Plane[3])+str(Plane[4])
        gc.collect()
        return S;
    
    def form_Ad(self):
        lane = []
        for z in range(0,int(self.w)):
            lane.append("")
        sheet = []
        for y in range(0,5):
            sheet.append(lane[:])
        Ad = []
        for x in range(0,5):
            Ad.append(sheet[:])
        gc.collect()
        return Ad;

    def theta(self,A):
        row = []
        Ad = self.form_Ad()
        for z in range(0,self.w):
            row.append(0)
        c = []
        d = []
        for x in range(0,5):
            c.append(row[:])
            d.append(row[:])
        C = c[:]
        D = d[:]
        for x in range(0,5):
            for z in range(0,self.w):
                C[x][z] = int(A[x][0][z])^int(A[x][1][z])^int(A[x][2][z])^int(A[x][3][z])^int(A[x][4][z])
        for x in range(0,5):
            for z in range(0,self.w):
                D[x][z] = str(int(C[(x-1)%5][z])^int(C[(x+1)%5][(z-1)%self.w]))
        for x in range(0,5):
            for y in range(0,5):
                for z in range(0,self.w):
                    Ad[x][y][z] = str(int(A[x][y][z])^int(D[x][z]))
        gc.collect()
        return Ad;

    def rho(self,A):
        Ad = self.form_Ad()
        for z in range(0,self.w):
            Ad[0][0][z] = A[0][0][z]
        x,y = 1,0
        for t in range(0,24):
            for z in range(0,self.w):
                int(1)
                Ad[x][y][z]=str(A[x][y][int((z-( ( (t+1)*(t+2) )/2)) %self.w)])
            x,y = y,((2*x)+(3*y))%5
        gc.collect()
        return Ad;

    def pi(self,A):
        Ad = self.form_Ad()
        for x in range(0,5):
            for y in range(0,5):
                #print("x={}, y={},".format((x+(3*y))%5,y))
                for z in range(0,self.w):
                    Ad[x][y][z] = str(A[(x+(3*y))%5][x][z])
        gc.collect()
        return Ad;

    def chi(self,A):
        Ad = self.form_Ad()
        #print(A)
        for x in range(0,5):
            for y in range(0,5):
                for z in range(0,self.w):
                    #print("t={} te={}".format(t,te))
                    Ad[x][y][z]=str(int(A[x][y][z])^((int(A[(x+1)%5][y][z])^1) and int(A[(x+2)%5][y][z])))
        gc.collect()
        return Ad;

    def rc(self,t):
        if t%255 == 0:
            return 1;
        R = list("10000000")
        for i in range(1,t%255):
            R = ["0"]+R
            R[0] = str(int(R[0])^int(R[8]))
            R[4] = str(int(R[4])^int(R[8]))
            R[5] = str(int(R[5])^int(R[8]))
            R[6] = str(int(R[6])^int(R[8]))
            R = R[0:8]
        gc.collect()
        return R[0];

    def l(self,A,ir): #theres a bug in here
        Ad = self.form_Ad()
        for x in range(0,5):
            for y in range(0,5):
                for z in range(0,self.w):
                    Ad[x][y][z] = A[x][y][z]
        RC = self.w*["0"]
        for j in range(0,self.fl+1):
            RC[(2**j)-1]=self.rc(j+(7*ir))
        for z in range(0,self.w):
            Ad[0][0][z]=str(int(Ad[0][0][z])^int(RC[z]))
        gc.collect()
        return Ad;

    def rnd(self,A,ir):
        #print(A[0][0])
        out = self.l(self.chi(self.pi(self.rho(self.theta(A)))),ir)
        return out;

    def KECCAKp(self,s): #s is string, self.nr = number of rounds self.b = length of s, w,fl are extra
        A = self.S_to_A(s)
        #print(A[0][0])
        for ir in range(12+(2*self.fl)-self.nr,12+(2*self.fl)-1):
            self.iter+=1
            #if self.iter>29998:
                #print("Iterations={},State={},b={},ir={}")
            A = self.rnd(A,ir)
        Sd = self.A_to_S(A)
        gc.collect()
        return Sd;

    def KECCAKf(self,s):
        self.nr = 12+(2*self.fl)
        gc.collect()
        return self.KECCAKp(s);
             
    def pad10n1(self,x,m):
        j = (-m-2)%x
        p=""
        for i in range(0,j):
            p=p+"0"
        p = "1"+p+"1"
        gc.collect()
        return p;

    def SPONGE(self,f,pad,r,N,d): #f = KECCAKp, pad = pad10n1, r = rate, N = string, d = int>0 = output length
        P=N+pad(r,len(N))
        n = len(p)/r
        c = self.b-r
        a = 0
        p = []
        for i in range(0,len(P),r):
            a=i
            p.append(P[i:i+r-1])
        del P
        del N
        del r
        gc.collect()
        S = ""
        for i in range(0,self.b):
            S=S+"0"
        e = ""
        gc.collect()
        for i in range(0,c):
              e=e+"0"
        for i in range(0,n-1):
            S = f(S^(p[i]+e))
        gc.collect()
        Z = ""
        while True:
            Z = Z+S[:r]
            if d <= len(Z):
                return Z[:d];
            gc.collect()
            S = f(S)

    def KECCAK(self,c,N,d): #c = capacity (will use 256),  N is string, d is length of output
        return SPONGE(1600,24,self.KECCAKp,self.pad10n1,1600-c,N,d);

class PRNG:
    def __init__(self,r,c,k,f,pa):
        self.r = r #rate (how many bit pullable)
        self.c = c #the number of bits not pullable
        self.k = k #sets the padding size
        self.f = f #function used during feed and fetch
        self.pa = pa #padding algorithm
        self.m = 0 #number of fetches since feed
        s = ""
        for i in range(0,self.r+self.c):
            s = s+"0"
        self.s = s #sets the state
        gc.collect()

    def main(self,request,inp):
        s = self.s
        f = self.f
        
        if request == "feed" and len(self.pad(inp))==self.k*self.r:
            print("FEED")
            p = []
            P = self.pad(inp)
            for i in range(0,len(P),self.r):
                p.append(P[i:i+self.r])
            ex = ""
            for i in range(0,self.c):
                ex = ex+"0"
            S = list(s)
            for i in range(0,self.k):
                gc.collect()
                P = p[i]
                use = P+ex
                for a in range(0,len(s)):
                    S[a]=str(int(s[a])^int(use[a]))
                s = "".join(S)
                #print("s-pre=",s)
                s = str(f(s))
                #
                print("s-post=",s)
            self.m = 0
            self.s = s
            gc.collect()
            return "";

        
        if request == "fetch":
            gc.collect()
            return self.squeeze(inp);

    def squeeze(self,l):
        # a is the number of available bits, s is the state, m is the number of bits gained, 
        s = self.s
        m = self.m
        r = self.r
        f = self.f
        if self.m == 0:
            a = self.r
        else:
            a = (0-self.m)%self.r
        full = ""
        t = 1
        print("start s =",s)
        gc.collect()
        while l>0:
            gc.collect()
            #print("squeze stage=",t);t+=1
            if a == 0:
                #print("s=",s)
                s = self.f(s)
                #print("s=",s)
                a = self.r
            ld = min(a,l)
            out = s[self.r-a:self.r-a+ld-1]
            a-=ld
            l-=ld
            m+=ld
            full=full+out
            self.s = s
        gc.collect()
        return full;

    def pad(self,inp):
        need = self.k+self.r - len(inp)
        x = len(inp)
        m = (x*need)-2
        a = self.pa(x,m)
        out = inp+a
        if len(out) != need:
            len_ = self.k*self.r
            out = out[:-1]+((len_-len(out))*"0")+"1"
            gc.collect()
        return out;

def runner(request,inp,sponge,gen):
    out = gen.main(request,inp)
    if out != "":
        return out;
    
def main(data_set,num_needed):
    k = 1000
    sponge = Sponge(1344+256)
    gen = PRNG(1344,256,k,sponge.KECCAKf,sponge.pad10n1)
    for i in data_set:
        gc.collect()
        runner("feed",i,sponge,gen)
    results = []
    print("###Generating###")
    for i in range(num_needed):
        gc.collect()
        out = ""
        for t in range(0,10):
            gc.collect()
            out = out+(runner("fetch",100000,sponge,gen))
        print(out)
        results.append(out)
    return results;
