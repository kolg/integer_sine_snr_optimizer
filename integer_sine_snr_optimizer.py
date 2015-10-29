#!/usr/bin/env python 

#title:		integer_sine_snr_optimizer.py
#description:	Optimizes the SNR for an integer sine with a given length by varying the amplitude
#author:	Martin Beha
#date:		20151029
#usage:		python integer_sine_snr_optimizer.py
#license:	GNU GPL v3.0
#website:	https://mbeha.wordpress.com/2015/10/30/snr-optimized-rounding-of-known-signals/
#==============================================================================

#import needed modules
import numpy as np
import matplotlib.pyplot as plt



#Variables that can be modified by the user are between the lines:
#------------------------------------------------

n = 24 			#sine with n values
Alimit = 2048 		#value that should not be reached
Amin = 0.5		#minimum value

rangeincrement = 0.25	#increment for the phase phi from zero (0) to half a value distance (1)

sinplotamplitude = 7	#amplitude for the sine to be plotted in the second figure
phiplot = 1		#phase to be plotted in the second figure


#------------------------------------------------

#variables needed in the following code
snr_db_maxarr=[]
snr_db_arr=[]
aopt=[]
asinopt=[]
amaxint=[]
snr_db_tempmax = 0
iterations = 0
iteration_phichange = [0]

Aopt_absolute = 0
Aopt_rounded_absolute = 0
phiopt_absolute = 0
snropt_absolute = -10

x = np.arange(n)
normsin = np.sin(2*np.pi*x/n)
epsilon = 1e-6	

#check rangeincrement
if rangeincrement > 1 or rangeincrement < 0:
  rangeincrement = 0.5
if rangeincrement == 0:
  rangeincrement = 1

#check if maximum phi should be 1 or 0.5 (values > 0.5 are redundant in some cases)
if (n % 2 > 0) and ((0.5*(1+epsilon)) % rangeincrement < epsilon):
  rangemax = 0.5+epsilon	#if n is odd and one value of phi is 0.5, all following values will be redundant
else:
  rangemax = 1+epsilon	#gerade
if rangeincrement == 1 or rangeincrement == 0:
  phiarr = [0]
else:
  phiarr = np.arange(0,rangemax,rangeincrement)

#loop to calculate all phases
for phi in phiarr:
  normsin = np.sin(2*np.pi*x/n + np.pi/n*phi)
  A = Amin
  if A < 0.5 or A >= Alimit:
    A = 0.5
    Amin = 0.5
  snr_db_tempmax = -10
  #loop to increase amplitude step by step until Alimit is reached
  
  #print status
  print "optimizing with phi=%1.4f" % (phi)
  
  while 1:
    #create rounded sine values
    sinval = A * np.sin(2*np.pi*x/n + np.pi/n*phi)
    sinvalround = np.ma.round(sinval)
    #check if limit is reached
    if max(sinvalround)>=Alimit:
      iteration_phichange.append(iterations)
      break
    if max(sinvalround)<1:
      A = A+0.01
      continue
    
    #get amplitude of the fundamental of the rounded sine
    Asin = np.sum(sinvalround*normsin)*2/n 
    
    #calculate residues and SNR
    rerror = Asin * normsin - sinvalround
    rerrorrms = np.sqrt(np.sum(np.power(rerror,2)))
    asinvalrms = np.sqrt(np.sum(np.power(Asin * normsin,2)))
    snr = asinvalrms / rerrorrms
    snr_db = 20*np.log10(snr)
    
    #memorize temporary maximum and other values for plotting
    if snr_db > snr_db_tempmax:
      snr_db_tempmax=snr_db
      if snr_db > snropt_absolute:
	Aopt_absolute = A
	Aopt_rounded_absolute = Asin
	phiopt_absolute = phi
	snropt_absolute = snr_db
    snr_db_maxarr.append(snr_db_tempmax)
    snr_db_arr.append(snr_db)
    aopt.append(A)
    asinopt.append(Asin)
    amaxint.append(max(sinvalround))
        
    iterations = iterations+1
    
    #calculate the next possible value of A to reach a different rounded sine
    nextstepdist = 1
    for i in range(n):
      if np.abs(normsin[i]) > 1e-6:
	if (sinvalround[i]+0.5*np.sign(sinval[i])-sinval[i])*A/sinval[i] < nextstepdist:
	  nextstepdist = (sinvalround[i]+0.5*np.sign(sinval[i])-sinval[i])*A/sinval[i]
    A = A + nextstepdist + epsilon
    
    
#print status and results
print "optimization finished\n"
print "Optimum value for N=%d in the range A=%1.1f...%1.1f:" % (n, Amin, Alimit)
print "phi=%1.1f, Aopt=%1.4f, Aopt_rounded=%1.4f, SNR_opt=%1.1fdB\n" % (phiopt_absolute, Aopt_absolute, Aopt_rounded_absolute, snropt_absolute)


#plot figure 1
plt.figure(figsize=(10,6)) 
plt.title('Signal to Quantization Noise Ratio over Amplitude')
plt.xlabel('Amplitude')
plt.ylabel('SNR / dB')

#color palette for values of phi
color=iter(plt.cm.rainbow(np.linspace(0,1,len(phiarr))))
for phinum in range(len(iteration_phichange)):
  if phinum == len(iteration_phichange)-1:
    break
  
  #get the corresponding values out of the array
  plotstart = iteration_phichange[phinum]
  plotend = iteration_phichange[phinum+1]

  col=next(color)
  
  #plot SNR
  plt.plot(asinopt[plotstart:plotend],snr_db_arr[plotstart:plotend],c=col,label="phi%d=%1.4f" % (phinum, phiarr[phinum]))
  
  #plot temporary maximum SNR
  plt.plot(asinopt[plotstart:plotend],snr_db_maxarr[plotstart:plotend],'--',c=col)

#plot rule of thumb as comparation
plt.plot(np.arange(Amin,Alimit),1.76+6.02*(1+np.log2(np.arange(Amin,Alimit))),'k', label='Noise Model')

plt.grid(b=True, which="both",color='0.3',linestyle='--')#u'major')


#select linear or logarithmic x axis depending on the range of A
if Alimit / Amin > 10:
  plt.xscale('log')
  
plt.legend(loc='upper left',fontsize=12)

#create and plot text
textstring = ( 'sine signal, N = ' + str(n) )
plt.annotate(textstring, xy=(0.95, 0.05), xycoords='axes fraction', bbox={'facecolor':'white', 'alpha':1, 'pad':12}, horizontalalignment='right', verticalalignment='bottom')



#plot figure 2
plt.figure(figsize=(10,6))

#get data for plotting nice curves
xdetail = np.arange(0,n,0.005)
sinplot = sinplotamplitude*np.sin(2*np.pi*x/n + np.pi/n*phiplot)
sinplotdetail = sinplotamplitude*np.sin(2*np.pi*xdetail/n + np.pi/n*phiplot)
sinplotrealamplitude = abs(np.fft.fft(np.round(sinplot))[1])/n*2
sinplotdetailreal = sinplotdetail / sinplotamplitude * sinplotrealamplitude
sinplotreal = sinplot /sinplotamplitude * sinplotrealamplitude
sinploterror = sinplotreal-np.round(sinplot)
sinplotsquarederror = sinploterror**2

#get data for the text
sinplotsquarederrorsum = sum(sinplotsquarederror)
sinplotrealsquaredsum = sum(sinplotreal**2)
snr_plot = 10 * np.log10(sinplotrealsquaredsum / sinplotsquarederrorsum)

plt.title('Sine Quantization')

#plot the things
plt.plot(xdetail,sinplotdetail,'g',label='Sine')
markerline, stemlines, baseline = plt.stem(x,sinplot,'g')
plt.setp(markerline, 'markerfacecolor', 'g', 'markeredgewidth', 1, 'markersize', 7)
plt.setp(stemlines, 'linewidth', 0)

plt.plot(xdetail,np.round(sinplotdetail),'b',label='Rounded Sine')
markerline, stemlines, baseline = plt.stem(x,np.round(sinplot),'b')
plt.setp(markerline, 'markerfacecolor', 'b', 'markeredgewidth', 1, 'markersize', 7)
plt.setp(stemlines, 'linewidth', 0)

plt.plot(xdetail,sinplotdetailreal,'cyan',label='Fundamental of Rounded Sine')

plt.plot(xdetail,sinplotdetailreal-np.round(sinplotdetail),'orange',label='Quantization Error')
markerline, stemlines, baseline = plt.stem(x,sinploterror,'orange')
plt.setp(markerline, 'markerfacecolor', 'orange', 'markeredgewidth', 1, 'markersize', 7)
plt.setp(stemlines, 'linewidth', 0)

plt.plot(xdetail,(sinplotdetailreal-np.round(sinplotdetail))**2,'magenta',label='Squared Error')
markerline, stemlines, baseline = plt.stem(x,sinplotsquarederror,'magenta')
plt.setp(markerline, 'markerfacecolor', 'magenta', 'markeredgewidth', 1, 'markersize', 7)
plt.setp(stemlines, 'linewidth', 0)

plt.setp(baseline, 'color','k', 'linewidth', 1)
plt.grid(b=True, which="both",color='0.3',linestyle='--')#u'major')
plt.xticks(np.arange(0, n, 1))

#fix y axis if needed for a nice figure
#plt.ylim(ymin = -4.5)
#plt.ylim(ymax = 4.5)

plt.legend(loc='upper right',fontsize=12)

#create and plot text
textstring = ( 'N = ' + str(n) \
	     + '\nphi = ' + ("%1.4f" % (phiplot)) \
	     + '\nA = ' + ("%1.2f" % (sinplotamplitude)) \
	     + '\nA_rounded = ' + ("%1.2f" % (sinplotrealamplitude)) \
	     + '\nSum of Squared Error: ' + ("%1.3f" % (sinplotsquarederrorsum)) \
	     + '\nSNR = ' + ("%1.1f" % (snr_plot)) + 'dB' )
plt.annotate(textstring, xy=(0.05, 0.05), xycoords='axes fraction', bbox={'facecolor':'white', 'alpha':1, 'pad':12})


#plot figure 3
plt.figure(figsize=(10,6)) 
plt.title('Signal to Quantization Noise Ratio over Maximum Integer Value')
plt.xlabel('Maximum Integer Value')
plt.ylabel('SNR / dB')

#color palette for values of phi
color=iter(plt.cm.rainbow(np.linspace(0,1,len(phiarr))))
for phinum in range(len(iteration_phichange)):
  if phinum == len(iteration_phichange)-1:
    break
  
  #get the corresponding values out of the array
  plotstart = iteration_phichange[phinum]
  plotend = iteration_phichange[phinum+1]

  col=next(color)
  
  #plot SNR
  plt.plot(amaxint[plotstart:plotend],snr_db_arr[plotstart:plotend],c=col,label="phi%d=%1.4f" % (phinum, phiarr[phinum]))
  
  #plot temporary maximum SNR
  plt.plot(amaxint[plotstart:plotend],snr_db_maxarr[plotstart:plotend],'--',c=col)
  
#plot rule of thumb as comparation
plt.plot(np.arange(Amin,Alimit),1.76+6.02*(1+np.log2(np.arange(Amin,Alimit))),'k', label='Noise Model')

plt.grid(b=True, which="both",color='0.3',linestyle='--')#u'major')

#select linear or logarithmic x axis depending on the range of A
if Alimit / Amin > 10:
  plt.xscale('log')
  
plt.legend(loc='upper left',fontsize=12)

#create and plot text
textstring = ( 'sine signal, N = ' + str(n) )
plt.annotate(textstring, xy=(0.95, 0.05), xycoords='axes fraction', bbox={'facecolor':'white', 'alpha':1, 'pad':12}, horizontalalignment='right', verticalalignment='bottom')

#show plots
plt.show() 
    
