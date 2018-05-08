from math import log

class Sponge:
    def __init__(self,b,nr=0):
        self.b = b
        self.w = self.b/25
        self.fl = int(log(self.w,2))
        if nr == 0:
            nr = 12+(2*self.fl)
        self.nr = nr
    def S_to_A(self,S):
        A  = form_Ad()
        for x in range(0,5):
            for y in range(0,5):
                for z in range(0,5):
                    A[x][y][z] = S[(w*((5*y)+x))+z]
    def A_to_S(self,A):
        d = []
        for a in range(0,5):
            d.append("")
        Lane = []
        for b in range(0,5):
            Lane.append(d)
        for i in range(0,5):
            for j in range(0,5):
                Lane[i][j] = A[i][j][0]+A[i][j][1]+A[i][j][2]+A[i][j][3]+A[i][j][4]
        Plane = ["","","","",""]
        for j in range(0,5):
            Plane[j]=Lane[0,j]+Lane[1,j]+Lane[2,j]+Lane[3,j]+Lane[4,j]
        S = Plane[0]+Plane[1]+Plane[2]+Plane[3]+Plane[4]
        return S;
    def form_Ad(self):
        lane = []
        for z in range(0,self.w):
            lane.append("")
        sheet = []
        for y in range(0,5):
            sheet.append(lane)
        Ad = []
        for x in range(0,5):
            Ad.append(sheet)
        return Ad;

    def theta(self,A):
        row = []
        Ad = form_Ad()
        for z in range(0,self.w):
            row.append(0)
        c = []
        d = []
        for x in range(0,5):
            c.append(row)
            d.append(row)
        for x in range(0,5):
            for z in range(0,self.w):
                C[x][z] = A[x][0][z]^A[x][1][z]^A[x][2][z]^A[x][3][z]^A[x][4][z]
        for x in range(0,5):
            for z in range(0,self.w):
                D[x][z] = C[(x-1)%5][z]^C[(x+1)%5][(z-1)%5]
        for x in range(0,5):
            for y in range(0,5):
                for z in range(0,self.w):
                    Ad[x][y][z] = A[x][y][z]^D[x][z]
        return Ad;

    def rho(self,A):
        Ad = form_Ad()
        for z in range(0,self.w):
            Ad[0][0][z] = A[0][0][z]
        x,y = 0,1
        for t in range(0,23):
            for z in range(0,self.w):
                Ad[x][y][z]=A[x][y][(z-((t+1)*(t+2))/2)%self.w]
                x,y = y,((2*x)+(3*y))%5
        return Ad;

    def pi(self,A):
        Ad = form_Ad()
        for x in range(0,5):
            for y in range(0,5):
                for z in range(0,self.w):
                    Ad[x][y][z] = A[(x+(3*y))%5][x][z]
        return Ad;

    def chi(self,A):
        Ad = form_Ad()
        for x in range(0,5):
            for y in range(0,5):
                for z in range(0,self.w):
                    Ad[x][y][z]=A[x][y][z]^(A[(x+1)%5][y][z]*A[(x+2)%5][y][z])
        return Ad;

    def rc(self,t):
        if t%255 == 0:
            return 1;
        R = "10000000"
        for i in range(1,t%255):
            R = "0"+R
            R[0] = str(int(R[0])^int(R[8]))
            R[4] = str(int(R[4])^int(R[8]))
            R[5] = str(int(R[5])^int(R[8]))
            R[6] = str(int(R[6])^int(R[8]))
            R = R[0:8]
        return R[0];

    def l(self,A,ir):
        Ad = form_Ad()
        for x in range(0,5):
            for y in range(0,5):
                for z in range(0,self.w):
                    Ad[x][y][z] = A[x][y][z]
        RC = ""
        for i in range(0,self.w):
            RC = RC+"0"
        for j in range(0,self.fl):
            RC[(2**j)-1]=rc(j+7*ir)
        for z in range(0,self.w):
            Ad[0][0][z]=int(Ad[0][0][0])^int(RC[z])
        return Ad;

    def rnd(self,A,ir):
        out = l(chi(pi(rho(theta(A)))),ir)
        return out;

    def KECCAKp(self,s): #s is string, self.nr = number of rounds self.b = length of s, w,fl are extra
        A = S_to_A(s)
        for ir in range(12+(2*self.fl)-self.nr,12+(2*self.fl)-1):
            A = rnd(A,ir)
        Sd = A_to_S(A)
        return Sd;

    def KECCAKf(self,s):
        self.nr = 12+(2*self.fl)
        return KECCAKp(s);
             
    def pad10n1(self,x,m):
        j = (-m-2)%x
        p=""
        for i in range(0,j):
            p=p+"0"
        p = "1"+p+"1"
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
        S = ""
        for i in range(0,self.b):
            S=S+"0"
        e = ""
        for i in range(0,c):
              e=e+"0"
        for i in range(0,n-1):
            S = f(S^(p[i]+e))
        Z = ""
        while True:
            Z = Z+S[:r]
            if d <= len(Z):
                return Z[:d];
            S = f(S)

    def KECCAK(self,c,N,d): #c = capacity (will use 256),  N is string, d is length of output
        return SPONGE(1600,24,KECCAKp,pad10n1,1600-c,N,d);

class PRNG:
    def __init__(self,r,c,k,f,pa):
        self.r = r
        self.c = c
        self.k = k
        self.f = f
        self.pa = pa
        self.m = 0
        s = ""
        for i in range(0,self.r+self.c):
            s = s+"0"
        self.s = s

    def main(self,request,inp):
        s = self.s
        f = self.f
        if request == "feed" and len(pad(inp))==self.k*self.r:
            p = []
            P = pad(inp)
            for i in range(0,len(P),self.r):
                p.append(P[i:i+self.r-1])
            ex = ""
            for i in range(0,self.c):
                ex = ex+"0"
            for i in range(0,self.k):
                P = p[i]
                use = P+ex
                for a in range(0,self.r+self.c):
                    s[a]=int(s[a])^int(use[a])
                s = f(s)
            self.m = 0
            self.s = s
            return "";
        if request == fetch:
            return self.squeeze(inp);

    def squeeze(self,l):
        s = self.s
        if self.m == 0:
            a = r
        else:
            a = (0-self.m)%self.r
        full = ""
        while l>0:
            if a == 0:
                s = self.f(s)
                a = self.r
            ld = min(a,l)
            out = s[self.r-a:self.r-a+ld-1]
            a-=ld
            l-=ld
            m+=ld
            full=full+out
        return full;

    def pad(self,inp):
        need = self.k+self.r - len(inp)
        x = len(inp)
        m = (x*need)-2
        """ (m+2)%x = need
             m+2 = x*need
             m = x*need - 2
        """
        a = self.pa(x,m)
        out = inp+a

def runner(request,inp,sponge,gen):
    out = gen.main(request,inp)
    if out != "":
        return out;
    
def main(data_set):
    k = 100
    sponge = Sponge(1344+256)
    gen = PRNG(1344,256,k,sponge.KECCAKf,sponge.pad10n1)
    for i in data_set:
        runner(feed,data,sponge,gen)
    out = runner(fetch,1000000,sponge,gen)
