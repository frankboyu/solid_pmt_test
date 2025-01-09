import os, sys
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

iteration = int(sys.argv[1])
group_list   = np.genfromtxt("list.txt", comments=None, dtype=str)
for group in group_list:
    run_list = np.genfromtxt("results/"+group+"/runs.txt",comments=None, dtype=str)[:,0]

    if(iteration==0):
        for run in run_list:
            directory = "benchtest/"+run+"/"
            for filename in os.listdir(directory):
                if filename.startswith("Fit"):
                    file_path = os.path.join(directory, filename)
                    os.remove(file_path)
    else:
        #NUMBER OF CHANNELS
        max_channel = 4

        #SAVE PLOTS AND PARAMETERS
        plots = []
        paras = np.zeros(24)

        for channel in range(1,max_channel+1):
            for run in run_list:
                this_plot = "benchtest/"+run+"/Fit_Ch0"+str(channel)+".png"
                plots.append(Image.open(this_plot))
                this_para = np.loadtxt("benchtest/"+run+"/Fit_Ch0"+str(channel)+".dat", comments=None, dtype=float)
                paras = np.vstack((paras, this_para))
        paras = paras[1:]

        plots[0].save("results/"+group+"/plots_iter0"+str(iteration)+"_10para.pdf", "PDF" ,resolution=100.0, save_all=True, append_images=plots[1:])
        np.savetxt("results/"+group+"/paras_iter0"+str(iteration)+"_10para.txt", paras, fmt='%i %i %i %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f')

        #AVERAGE THE SPE PARAMETERS IF IT'S THE FIRST ITERATION IN THE CYCLE
        if(iteration%2==1):
            for channel in range(1,max_channel+1):
                paras_new = paras[len(run_list)*(channel-1):len(run_list)*channel]

                paras_new[:,8]  = np.average(paras_new[:,8])
                paras_new[:,9]  = np.average(paras_new[:,9])
                paras_new[:,10] = np.average(paras_new[:,10])
                paras_new[:,11] = np.average(paras_new[:,11])
                paras_new[:,12] = np.average(paras_new[:,12])
                paras_new[:,13] = np.average(paras_new[:,13])

                for i, run in enumerate(run_list):
                    paras_print = paras_new[i,:].reshape(1, paras_new.shape[1])
                    np.savetxt("benchtest/"+run+"/Fit_Ch0"+str(channel)+".dat", paras_print, fmt='%i %i %i %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f')
