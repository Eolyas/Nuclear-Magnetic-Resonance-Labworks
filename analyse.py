import re
import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os


def exponential_fit(x,a,b):
    return a*np.exp(-x/b)

#Call the function with the folder name holding the .dat files. The files are assumed sorted in increasing delays.
def analysis(folder_name):
    filelist = []
    freq_dict = dict()
    for filename in os.listdir(folder_name):
        if re.findall(r'\d+.dat',filename):
            filelist.append(filename)
    #The delay is in power of 2 (2,4,8,16).
    #To modify the delay, comment the next line and uncomment the 24th line below.
    timelist = [2^(n+1) for n in np.size(filelist)]
    for i in range(np.size(filelist)):
        #timelist.append(input("delay for "+str(i+1)+"th file : "))
        file = open(folder_name+"/"+filelist[i],'r')
        lines = file.readlines()
        for j in range(np.size(lines)):
            temp = re.findall(r'\d+.\d+\s+-?\d+.\d+',lines[j])
            if temp!=[]:
                temp2 = re.split(r'\s+',temp[0])
                freq = temp2[0]
                amp = temp2[1]
                if not(freq in freq_dict):
                    freq_dict[freq] = []
                freq_dict[freq].append(amp)
    T2_list = []
    freq_list = []
    for i in freq_dict:
        params, covariance = curve_fit(exponential_fit,timelist,[float(x) for x in freq_dict[i]])
        fitted_a,T2 = params
        T2_list.append(float(T2))
        freq_list.append(float(i))
    fig,ax = plt.subplots()
    plt.plot(freq_list,T2_list)
    ax.set_xlabel("frequency")
    ax.set_ylabel("T2")
    plt.show() 
    return 0

def T2_analysis_deprecated():
    #file_list = ['24121942.dat','24121944.dat','24121945.dat','24121946.dat','24121947.dat']
    file_list = ['24121954.dat','24121952.dat','24121950.dat','24121948.dat','24121955.dat','24121956.dat']
    #time_list = [4,8,16,32,64]
    time_list = [0.2,0.3,0.6,1,1.3,1.6]
    integ = []
    nb_files = np.size(file_list)
    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))

    for j in range(nb_files):
        file = open("cofeet4copie/"+file_list[j], 'r')
        list_lines = file.readlines()
        size = np.size(list_lines)
        x = np.zeros(size)
        y = np.zeros(size)
        amp = np.zeros(size)
        for i in range(size):
            result = re.findall(r'-?\d+.\d+',list_lines[i])
            x[i] = result[0]
            y[i] = result[1]
            amp[i] = result[2]
        y*=np.cos(amp)
        y_normalized = (y - np.min(y))/(np.max(y)-np.min(y))
        ax1.plot(x,y,label="T="+str(time_list[j]) + "ms")
        ax1.set_xlabel('frequency')
        ax1.set_ylabel('amplitude')
        if time_list[j] == 256:
            var = round(np.var(y_normalized),3)
            ax2.plot(x,y_normalized,label="T="+str(time_list[j]) + "ms var = " + str(var),alpha=0.5)
        elif time_list[j] == 128:
            var = round(np.var(y_normalized),3)
            ax2.plot(x,y_normalized,label="T="+str(time_list[j]) + "ms var = " + str(var),alpha=0.8)
        else:
            var = round(np.var(y_normalized),3)
            ax2.plot(x,y_normalized,label="T="+str(time_list[j]) + "ms var = " + str(var))
        ax2.set_xlabel('frequency')
        ax2.set_ylabel('amplitude')
        integ.append(integrate.trapezoid(y,x))
        file.close
    ax1.legend()
    ax2.legend()
    integ = (integ-np.min(integ))/(np.max(integ)-np.min(integ))
    params, covariance = curve_fit(exponential_fit,time_list,integ)
    fitted_a,fitted_b = params
    fitted_y = np.zeros(nb_files)
    for i in range(nb_files):
        fitted_y[i] = exponential_fit(time_list[i],fitted_a,fitted_b)
    ax3.plot(time_list,integ,label="experimental curve")
    ax3.plot(time_list,fitted_y,label="fitted curve T2 = " + str(round(fitted_b,3)))
    ax3.legend()
    ax3.set_ylim(ymin = 0)
    ax3.set_xlabel('delay')
    ax3.set_ylabel('integrated frequencies')
    plt.savefig("graph T0 ambient")
    plt.show()
    return 0


def plot3D(folder_name):
    folder_name = str(folder_name)
    sin_file = ""
    for filename in os.listdir(folder_name):
        if re.findall(r'.sin',filename):
            sin_file = filename
    main_file = open(folder_name + "/" + sin_file,'r')
    main_lines = main_file.readlines()
    files_list = []
    field_list = []
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for i in range(np.size(main_lines)):
        if not i<6:
            temp = re.findall(r'\d+.dat\s+\d+',main_lines[i])
            temp2 = re.split(r'\s+',temp[0])
            files_list.append(temp2[0])
            field_list.append(float(temp2[1]))
    main_file.close()
    #sorting files and fields
    files_list = np.array(files_list)
    field_list = np.array(field_list)
    sorting_indices = np.argsort(field_list)
    sorted_fields = field_list[sorting_indices]
    sorted_files = files_list[sorting_indices]
    #
    amp_list = np.zeros((np.size(files_list),256))
    freq_list = np.zeros(256)
    for i in range(np.size(files_list)):
        file = open(folder_name + "/" + str(sorted_files[i]),'r')
        lines = file.readlines()
        k=0
        for j in lines:
            temp3 = re.findall(r'\d+.\d+\s+\d+.\d+',j)
            if temp3:
                temp4 = re.split(r'\s+',temp3[0])
                freq_list[k] = float(temp4[0])
                amp_list[i][k] = float(temp4[1])
                k+=1
    x,y = np.meshgrid(freq_list,[np.log(x) for x in sorted_fields])
    ax.plot_surface(x,y,amp_list,cmap="viridis")
    plt.show()
        
    
    return 0

plot3D("patapouf")