import CommonFunctions
import numpy as np
import matplotlib.pyplot as plt  

Name = "Buchnera_aphidicola"
plt.figure(figsize=(12, 14))  
ax = plt.subplot(311)  
plt.ylabel("No of substrings", fontsize=16)  
plt.xlabel("Read Length", fontsize=16) 
plt.xlim(0, 30)  
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
plt.semilogy(XAxis,YAxis,marker='s', lw=0.5, color="black", alpha=0.3)

plt.subplot(312)  
plt.ylabel("No of repeats", fontsize=16)  
plt.xlabel("Read Length", fontsize=16) 
plt.xlim(0, 30)  
plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on") 
Array = CommonFunctions.FiletoArray('Read_lengths'+Name+'.csv')
XAxis = []
for i in Array:
    XAxis += [ int(i[0]) ]
YAxis = []
for i in Array:
    YAxis +=  [int(i[2]) ]
    
plt.semilogy(XAxis,YAxis, marker='s', lw=0.5, color="blue", alpha=0.3)

plt.subplot(313)  
plt.ylabel("Uncertainty", fontsize=16)  
plt.xlabel("Read Length", fontsize=16) 
plt.xlim(0, 30)  
Array = CommonFunctions.FiletoArray('U_Summary_'+Name+'.csv')
XAxis = []
for i in Array:
    XAxis += [ int(i[0]) ]
YAxis = []
for i in Array:
    YAxis +=  [float(i[1]) ]
    
plt.semilogy(XAxis,YAxis, marker='s', lw=0.5, color="blue", alpha=0.3)

plt.show()
