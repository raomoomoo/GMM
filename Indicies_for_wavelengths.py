import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
import paramiko
from scp import SCPClient
import pwinput
import os
import sys
import getpass

#This file will actually craeate all the different wavelengths for set positions, and will be used to create the specific monomer, central only, and varies size for configuration of particles
#use user promopt for wavs, and aggs 

data = pd.read_csv('draine_optics.dat', header=None, delimiter=r"\s+", engine='python')
data.columns = ['Wav', 'Real', 'Img']

#data
#a nice plot of he refractive index 
# x=data['Wav']
# y=data['Img']
# plt.figure(figsize=(10,10))
# plt.rcParams.update({'font.size': 20})
# plt.plot(x,y, '.')
# plt.xlabel('wav')
# plt.ylabel('Img')
# plt.show()


wavs = np.arange(550, 2500, 50).tolist()

index = []
for i in wavs:
    f  = interpolate.interp1d(np.log(data.Wav), np.log(data.Img), fill_value='extrapolate')
    g  = interpolate.interp1d(np.log(data.Wav), np.log(data.Real), fill_value='extrapolate')
    im = np.exp(f(np.log(i)))
    rel= np.exp(g(np.log(i)))
    index.append([i,rel,im])

index
index_data=np.array(index)

#plot range of vale in interpolate and expolate in range 
#make sure is straight line
#check for ressonance 

# x=str(index_data[j, 0]/1000),
# y=index_data[:, 0]
# #plt.figure(figsize=(10,10))
# plt.rcParams.update({'font.size': 20})
# plt.plot(np.log(x),np.log(y), '.')
# plt.xlabel('real')
# plt.ylabel('wav')
# plt.show()

#A  function to edit the text file, and add stuff at the top  
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


data_Agg = pd.read_csv('/Users/raomorusupalli/Documents/UNI/Honours/project/ozstar_out_models/aggregate_260_3.k', header=None,skiprows=2,delimiter=r"\s+", engine='python' )
data_Agg.columns = ['x', 'y', 'z','size','Re','Img']

#just to look at data 
#data_Agg

# A function to check the progress of upload to remote machine 
def progress(filename, size, sent):
    sys.stdout.write("%s\'progress: %.2f%%   \r" % (filename, float(sent)/float(size)*100) )

#automating the transfer to a remote machine (oooh fancy !!)
ssh = paramiko.SSHClient()
ssh.load_host_keys(os.path.expanduser('~/.ssh/known_hosts'))

Pass = getpass.getpass('Enter your password: ')

ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('ozstar.swin.edu.au', username='rmorusup', password=Pass)
scp = SCPClient(ssh.get_transport(), progress = progress)   

#edits text file for the range of wavelengths from a precalculates aggregate model 
for i in wavs:
    for j in range(len(wavs)):
        data_Agg['Img']=index_data[j,2]
        data_Agg['Re']=index_data[j,1]
        filename = 'Aggregate_260_'+str(i)+'.k'
        np.savetxt(filename, data_Agg.values, fmt='%f' )
        list_of_lines = [str(261), str(i/1000)]
        [line_prepender(filename, line) for line in list_of_lines ]
        file_remote = r'/home/rmorusup/runs/Agg_models/AGG_260_3_wav/'+str(i)+'/'
        scp.put(filename, file_remote)
        
scp.close()
ssh.close()
           

#same as above for just 1 particle, a single size grain
#look into making this faster 
for i in wavs:
    for key in Agg:
        data_Agg['sizes'] = Agg[key]
        filename = str(key)+'_'+str(i)+'.k'
        for j in range(len(wavs)):
            data_Agg['Img']=index_data[j,2]
            data_Agg['Re']=index_data[j,1]
            np.savetxt(filename, data_Agg.iloc[0],newline=' ', delimiter=',', fmt='%f' )
            list_of_lines = [str(1), str(i/1000)]
            [line_prepender(filename, line) for line in list_of_lines]
            file_remote = r'/home/rmorusup/runs/Agg_models/AGG_260_3_wav/'+str(i)+'/'
            scp.put(filename, file_remote)
        
scp.close()
ssh.close()