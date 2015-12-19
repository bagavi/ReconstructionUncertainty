import CommonFunctions
import numpy as np
import matplotlib.pyplot as plt  

Name = "RhodobacterSphaeroides"
x_start = 0
x_end = 700
y_start = 0
y_end = 2^22
    
plt.figure(figsize=(12, 14))  
ax = plt.subplot(311)  
plt.ylabel("No of substrings", fontsize=16)  
plt.xlabel("Read Length", fontsize=16) 
plt.xlim(x_start, x_end)  
#plt.ylim(y_start, y_end)  
#plt.yticks(np.arange(y_start,y_end,2))

plt.title(Name)
ax.get_xaxis().tick_bottom()    
ax.get_yaxis().tick_left()    
Array = CommonFunctions.FiletoArray('Read_lengths'+Name+'.csv')
XAxis = []
for i in Array:
    XAxis += [ int(i[0]) ]
YAxis = []
for i in Array:
    YAxis +=  [int(i[1]) ]
print(len(XAxis),len(YAxis))    
plt.semilogy(XAxis,YAxis,marker='s', lw=0.5, color="black", basey=2, alpha=1)

plt.subplot(312)  
plt.ylabel("No of repeats", fontsize=16)  
plt.xlabel("Read Length", fontsize=16) 
plt.xlim(x_start, x_end)  
#plt.ylim(y_start, y_end)  

plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on") 
Array = CommonFunctions.FiletoArray('Read_lengths'+Name+'.csv')
XAxis = []
for i in Array:
    XAxis += [ int(i[0]) ]
YAxis = []
for i in Array:
    YAxis +=  [int(i[2]) ]
    
plt.semilogy(XAxis,YAxis, marker='s', lw=0.5, color="blue", basey=2,alpha=1)
 
plt.subplot(313)  
plt.ylabel("Uncertainty", fontsize=16)  
plt.xlabel("Read Length", fontsize=16) 
plt.xlim(x_start, x_end) 
#plt.ylim(y_start, y_end)  
 
Array = CommonFunctions.FiletoArray('U_Summary_'+Name+'.csv')
XAxis = []
for i in Array:
    XAxis += [ int(i[0]) ]
YAxis = []
for i in Array:
    YAxis +=  [float(i[1]) ]
     
plt.semilogy(XAxis,YAxis, marker='s', lw=0.5, color="blue", basey=2,alpha=0.3)
 
YAxis = []
for i in Array:
    if float(i[2]) < 1:
        YAxis += [0]
    else:    
        YAxis +=  [float(i[2]) ]
 
plt.semilogy(XAxis,YAxis, marker='s', lw=0.5, color="red", basey=2,alpha=1)


plt.show()
