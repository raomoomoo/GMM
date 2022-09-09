import random as ran
import numpy as np
from math import sqrt
import cProfile
import time
from functools import cache
import os.path
import sys
import gc


n_m = int(sys.argv[1])

r_mm = float(sys.argv[2])
r_m = r_mm/1e3
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

# while True: #choose monomer size
#     try:
#         r_m = float(input('Enter size of monomers in mm e.g 0.01: '))
#     except ValueError:
#         # Not a valid number
#         print('You must enter a vaild number')
#         r_m = float(input('Enter the number of monomers: '))
#     else:
#         # No error; stop the loop
#         break

#setup arrays and initialise values 

x_a, y_a, z_a = [],[],[]

Wavelength = 0.870
Re = 3.408100
Im = 5.6983002E-02

#n_m =   number of monomeres
r_c = 1 
#r_m  = 0.01

fac = (r_c + r_m) * (1+1e-6)

save_path = '/Users/raomorusupalli/Documents/UNI/Honours/project/GMM/'

name_of_file ='aggregate_'+str(n_m)

completeName = os.path.join(save_path, name_of_file+".k")

start_time = time.time()
np.seterr(invalid='ignore')


@cache
def main():
    count = 0
    tried =0
    actp = []
    t1 = time.time()
    with open(completeName,'w') as f:
        f.write(str(Wavelength)+'\n')
        f.write(str(n_m+1)+'\n')
        f.write(" ".join([str(0.), str(0.), str(0.), str(r_c), str(Re), str(Im), '\n'])) #initial position of centre monomer
        while count < n_m: # i and j are indicies for 2 points to compare, j starts at i + 1 because all values of j=0, ..., i will be compared are a different iteration
                

            pos = ran.uniform(0.,1.)
            z = 2.0 * pos - 1.0
            r = np.sqrt(1-z*z)

            pos = ran.uniform(0.,1.)
            ARGMT = np.pi*(2.0*pos- 1.0 )
            x = r * np.cos(ARGMT)
            y = r * np.sin(ARGMT)

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
       
                    # d = np.sqrt(np.dot(posvec,actp[-1]))
                    # if fac/d > 1:
                    #     overlap = True
                    #     break 
                
        
            if not overlap:
                gc.disable()    
                actp.append(posvec)

                #print (posvec)
                count+=1
               # t1 = time.time()
                print(count)   
            else:
                tried+=1
       
        
        for val in actp:
           f.write('{} {} {}'.format(val[0], val[1], val[2]) +' '+str(r_m)+' '+str(Re)+' '+str(Im)+'\n')            
        f.close()

if __name__ == '__main__':

     cProfile.run('main()')    
     print("--- %s seconds ---" % (time.time() - start_time))