import sys
import math
import numpy as np
import symnmf as sp
from sklearn.metrics import silhouette_score

np.random.seed(0)
iter=300
epsilon=0.0001

def main() :
                       
    #KMEANS:
    centroids=[0 for i in range(k)]
    itnum=0
    y_1=[0 for i in range(n)]
   
    for i in range(k):
        centroids[i]=points[i]
    
   
    
    zero=[0 for i in range(len(points[0]))]
    diff=[]
    for i in range(k):
        diff.append(distance(zero,centroids[i]))       

    
    while(max(diff)>=0.0001 or itnum<=iter):
 
       clusters=[[] for i in range(k)]
       
       #cluster building
       for i in range(len(points)):
           mincentroid= min(centroids,key=lambda e: distance(e,points[i]))
           minindex=centroids.index(mincentroid)
    
           clusters[minindex].append(points[i])
           y_1[i]=minindex   
        #centroid updating
       centroidsprev=centroids.copy()
           
       for i in range(k):
           a = len(points[0])
           sumlst=zero.copy()
           
           for point in clusters[i]:
               sumlst=sumpoints(sumlst,point)
                 
           centroids[i]=dividepoint(sumlst,len(clusters[i])) 
          
           
        #centroid difference 
        
       for i in range(k):
            diff[i]= distance(centroidsprev[i],centroids[i])
        
       itnum+=1
    
    y_1=np.array(y_1)
    p=np.array(points)
    r_1=silhouette_score(p,y_1)
    
    #SYMNMF:
    W=sp.norm(points)
    W_np=np.array(W)
    m=np.mean(W_np)
    H=np.random.uniform(0,2*math.sqrt(m/k),(n,k))
    H=H.tolist()
    result=sp.symnmf(H,W,k)
    
    y_2=[0 for i in range(n)]
    
    for i in range(len(result)):
        y_2[i]=result[i].index(max(result[i]))
        
    y_2=np.array(y_2)
    r_2=silhouette_score(p,y_2)
    
    #results:
    print("nmf:","{:.4f}".format(r_2))
    print("kmeans:","{:.4f}".format(r_1))
          
#Eclidean distance    

def distance(a,b):
    sum=0

    for i in range(len(a)):
        sum= sum+ pow((a[i] - b[i]),2)
    
    return math.sqrt(sum)   

def sumpoints(a,b):

    c=[0 for i in range(len(a))]
    for i in range(len(a)):
        c[i]=a[i]+ b[i]
        
    return c 

def dividepoint (a,c):
    b=[]


    for i in range(len(a)):
        b.append(float(a[i]/c))
    return b


#Arguments Reading:
args=sys.argv[1:]
k=int(args[0])
file=open(args[1],"r")

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