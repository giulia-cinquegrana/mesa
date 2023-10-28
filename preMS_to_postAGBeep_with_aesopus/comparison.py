import matplotlib.pyplot as plt
import numpy as np
import mesa_reader as mr

aesopus = mr.MesaData("LOGS/topostAGBeep/history.data")
freedman = mr.MesaData("../preMS_to_postAGBeep_with_aesopus_freedman/LOGS/topostAGBeep/history.data")

fig=plt.figure()

ax=fig.add_subplot(111,label='1')
ax2=fig.add_subplot(111, label='2', frame_on=False)

ax.plot(aesopus.star_age, aesopus.log_LHe, color='#332288', label="aesopus")

ax.axis(xmin=780197000,xmax=781357000) 
ax.xaxis.set_label_position('bottom') 
ax.xaxis.tick_bottom()

ax.yaxis.set_label_position('left') 
ax.yaxis.tick_left()

ax.tick_params(axis='x', colors='#332288')

#ax.set_xlabel('Stellar age [yrs]', color='#004488', fontsize=15)
#ax.set_ylabel(r'Log(L$_{\rm He}$[L$_{\odot}$])', fontsize=15)

#------------------------------------------------------------------------------

ax2.plot(freedman.star_age, freedman.log_LHe, color='#BB5566', label="freedman")

ax2.axis(xmin=781240000,xmax=7.824e8) 
ax2.xaxis.set_label_position('top') 
ax2.xaxis.tick_top()

ax2.yaxis.set_label_position('right') 
ax2.yaxis.tick_right()

ax2.tick_params(axis='x', colors='#BB5566')
ax.legend()

plt.show()

#ax.set_xlabel('Stellar age [yrs]', color='#004488', fontsize=15)
#ax.set_ylabel(r'Log(L$_{\rm He}$[L$_{\odot}$])', fontsize=15)