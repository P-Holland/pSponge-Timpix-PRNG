"""Will use the first and last 5 files from lucid to use a inital seed
Will also collect pick random files from the remaining 80 file batches taking 5 from each
Total of 410 files as seed, feeding 2 every 5 fetches
 - 10*feed ,5* fetch, feed , ...
"""
import struct
import Generator
#import SHA3
from Tests import tests #main is to run
import os
def binary(num):
    return bin(struct.unpack('!i',struct.pack('!f',1.0))[0])
def fetch_file(num): #gets the file the number num and gets the usefull part I need
    location = os.path.dirname(__file__)
    file = open((location+"/files/lucid{}.txt".format(num)),"r")
    a = []
    for i in range(0,256**2):
        a.append("0")
    b = file.readlines()
    for i in b:
        temp = i.split("\t")
        if float(temp[2])>1:
            t = float(temp[2].strip("\n"))
            t = int(round(t))
            a[(int(temp[0])*256)+int(temp[1])]=(bin(t)[2:])
        else:
            a[(int(temp[0])*256)+int(temp[1])] = (binary(temp[2].strip("\n"))[2:])
    last = 2
    use = ""
    for i in a:
        if i!=last:
            last = i[0]
            use=use+i
    return use;

def main():
    fails=0
    for i in range(8,10):
        full = fetch_file(i)
        use = full[:1344*1000]
        print(use)
        nums=Generator.main([use],10)
            #nums = []
            #for t in range(10):
             #   nums.append(SHA3.SHAKE256([int(use)],10**6))
        count = 0
        file=open("nums{}.txt".format(i),"a")
        for a in nums:
            file.write(a)
        file.close()
        print("TESTING")
        for ep in nums:
            count+=1
            print("Testing number {}".format(count))
            test = tests(ep,0.01)
            result = test.main()
            if not result:
                print("!!!file {} number {} Failed!!!".format(i,count))
            else:
                print("!!!file {} number {} Passed!!!".format(i,count))
main()
                
