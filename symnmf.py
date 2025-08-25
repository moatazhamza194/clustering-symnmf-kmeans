import sys
import math
import numpy as np
import symnmf as sp

np.random.seed(0)
max_iter=300
epsilon=0.0001

def main():
    result=[[]]
    if(goal=="sym"):
        result=sp.sym(points)
    
    
    if(goal=="ddg"):
        result=sp.ddg(points)
    
        
    if(goal=="norm"):
        result=sp.norm(points)
        
    
    if(goal=="symnmf"):
        #Building initial H:
        W=sp.norm(points)
        W_np=np.array(W)
        m=np.mean(W_np)
        H=np.random.uniform(0,2*math.sqrt(m/k),(n,k))
        H=H.tolist()
        result=sp.symnmf(H,W,k)
        
    for i in range(len(result)):
        print(','.join('{0:.4f}'.format(x) for x in result[i]))    

        
    return 0

#Arguments Reading:
args=sys.argv[1:]
k=int(args[0])
goal=args[1]
file=open(args[2],"r")

points=[]

#File Reading:
curr_line = 0
for line in file:
    str_line = line.split(",")
    int_line = []
    for i in range(len(str_line)):
        int_line.append(float(str_line[i]))
    points.append(int_line)  # lines list[i] will contain the i'th vector in the list
    if line != "\n":
        curr_line += 1
        
n=len(points)
dim=len(points[0])
main()       



