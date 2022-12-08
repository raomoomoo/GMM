import random as ran
import numpy as np
from math import sqrt
import cProfile
import time
from functools import cache
from itertools import chain
import os.path
import sys
import gc

# while True: #choose monomer number 
#     try:
#         n_m = int(input('Enter the number of monomers: '))
#     except ValueError:
#         # Not a valid number
#         print('You must enter a vaild number')
#         n_m = int(input('Enter the number of monomers: '))
#     else:
#         # No error; stop the loop
#         break


while True: #choose monomer size
    try:
        r_mm = float(input('Enter the size of monomers: '))
    except ValueError:
        # Not a valid number
        print('You must enter a vaild number')
        r_mm = float(input('Enter the number of monomers: '))
    else:
        # No error; stop the loop
        break

while True: #choose layers
    try:
        l = int(input('number of layers: '))
    except ValueError:
        # Not a valid number
        print('You must enter a vaild number')
        l = int(input('number of layers: '))
    else:
        # No error; stop the loop
        break

        
#create a list of empty lists of size potentially do this part later after looping through the size of the arrays
#, specify dtype=object
#np.empty([10, 4], dtype=object)

layers = np.empty(l, dtype=object)
#print(layers)

# loop through a an input satement to determine the length of an empty list
#where the empty list is determined by l 
for i in range(l):
    #crea
    n=int(input('how many particles in layer:?'))
    layers[i] =[0] * n
#^figure out how to wrap the above in so you know what layer it is 


#print(layers)

n_m= int(sum(len(x) for x in layers))
#chamge 'total' = n_m
#print(total)


#break 
x_a, y_a, z_a = [],[],[]

Wavelength = 0.870
Re = 3.408100
Im = 5.6983002E-02

#n_m =   number of monomeres
r_c = 1 
r_m = r_mm/1e3

#fac = (r_c + r_m) * (1+1e-6)

save_path = '/Users/raomorusupalli/Documents/UNI/Honours/project/GMM/'

name_of_file ='aggregate_'+str(n_m)+'_'+str(l)

completeName = os.path.join(save_path, name_of_file+".k")

start_time = time.time()
np.seterr(invalid='ignore')


def main():
    count = 0
    tried =0
    j=0 # counter for the layer 
    fill=0 #the amount of particles in the given layer 
    actp = []
    t1 = time.time()
    with open(completeName,'w') as f:
        f.write(str(Wavelength)+'\n')
        f.write(str(n_m+1)+'_'+str(l)+'\n')
        f.write(" ".join([str(0.), str(0.), str(0.), str(r_c), str(Re), str(Im), '\n'])) #initial position of centre monomer
       # while count < n_m: # i and j are indicies for 2 points to compare, j starts at i + 1 because all values of j=0, ..., i will be compared are a different iteration 
        for i in range(l):
            for j in layers[i]:
                while fill <= len(layers[i])-1: #put the particles in the given layer  


                    pos = ran.uniform(0.,1.)
                    z = 2.0 * pos - 1.0
                    r = np.sqrt(1-z*z)

                    pos = ran.uniform(0.,1.)
                    ARGMT = np.pi*(2.0*pos- 1.0 )
                    x = r * np.cos(ARGMT)
                    y = r * np.sin(ARGMT)

                    #This is the stuff the actualy adds particles to the surface via the position vector, add multiple layer by changing 
                    #this something like fac = (r_c + l*r_m)
                    #where l is layer in the code 
                    fac = (r_c + (i+1)*r_m)
                    x *= fac
                    y *= fac
                    z *= fac

                    posvec = np.array([x,y,z])
                    overlap = False


                    if not overlap:
                          for p in actp:
                            # remove the doubling of the posvec[.]-p[.] in the below statement
                            # find a a way to make the staements absolute ie a == b  instead of a < b.
                            gc.disable()   
                            xt = np.abs(posvec[0]-p[0]) - 2*r_m
                            yt =  np.abs(posvec[1]-p[1]) - 2*r_m
                            zt = np.abs(posvec[2]-p[2]) -  2*r_m
                            rt = (posvec[0]-p[0])**2+(posvec[1]-p[1])**2+(posvec[2]-p[2])**2 - 4*r_m**2
                            stat = [abs(xt) != xt, abs(yt) != yt,abs(zt) != zt,abs(rt) != rt]
                            if all(stat):
                                overlap = True
                                break


                    if not overlap:
                        gc.disable()    
                        actp.append(posvec)
                        fill+=1
                        #print (posvec)
                       # count+=1
                       # t1 = time.time()
                        print(fill)   
                    else:
                        tried+=1
            if fill > len(layers[i])-1:
                fill=0
                j+=1
        
        for val in actp:
            f.write('{} {} {}'.format(val[0], val[1], val[2]) +' '+str(r_m)+' '+str(Re)+' '+str(Im)+'\n')            
        f.close()

if __name__ == '__main__':
    main()    
    print("--- %s seconds ---" % (time.time() - start_time))