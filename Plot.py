import CommonFunctions
import numpy as np
import matplotlib.pyplot as plt  

plt.figure(figsize=(12, 14))  
plt.subplot(311)  
plt.ylabel("No of substrings", fontsize=16)  
plt.xlabel("Read Length", fontsize=16) 
plt.title("Buchnera_aphidicola")
Array = CommonFunctions.FiletoArray('Read_lengthsBuchnera_aphidicola.csv', Int=True, Offset = 1)
Array = np.array(Array) 
plt.semilogy(Array[:,[0]], Array[:,[1]], "--", lw=0.5, color="black", alpha=0.3)

plt.subplot(312)  
plt.ylabel("No of repeats", fontsize=16)  
plt.xlabel("Read Length", fontsize=16) 
plt.title("Buchnera_aphidicola")
Array = CommonFunctions.FiletoArray('Read_lengthsBuchnera_aphidicola.csv', Int=True, Offset = 1)
Array = np.array(Array) 
plt.semilogy(Array[:,[0]], Array[:,[2]], "--", lw=0.5, color="black", alpha=0.3)
# 
# plt.subplot(313)  
# plt.ylabel("No of substrings", fontsize=16)  
# plt.xlabel("Read Length", fontsize=16) 
# plt.title("Buchnera_aphidicola")
# Array = CommonFunctions.FiletoArray('Read_lengthsBuchnera_aphidicola.csv', Int=True, Offset = 1)
# Array = np.array(Array) 
# plt.semilogy(Array[:,[0]], Array[:,[1]], "--", lw=0.5, color="black", alpha=0.3)



plt.show()
