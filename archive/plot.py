import numpy as np
import matplotlib.pyplot as plt
from os import path


files = [    
['20230802_1942', '100%'],
['20230802_1943', '80%'],
['20230802_1944', '60%'],
['20230802_1945', '40%'],
['20230802_1946', '20%'],
]


plt.figure(figsize=(20,10))
for (i, [name, tag]) in enumerate(files):

    if (0<=i<=4):
        color='r'
    elif (5<=i<=9):
        color='g'
    elif (10<=i<=14):
        color='b'
    else:
        color='m'        
    
    for j in range(1,5):
        data = np.loadtxt("data/"+name+"/Histo_Ch0"+str(j)+".txt") 
        label = " Ch" + str(j) +"," + tag
        plt.plot(data[:,0], data[:,1]/data[:,1].sum(), label=label)
#         plt.plot(data[:,0], data[:,1]/data[:,1].sum(), label=label, c=color)
#         plt.plot(data[:,0], data[:,1], label=label, c=color)
        plt.yscale('log')
        plt.xlim(100,1000)
        plt.legend(fontsize=12)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        print(name)
plt.savefig('foo.png',bbox_inches='tight')
