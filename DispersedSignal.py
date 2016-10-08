import numpy as np
import sys 
import matplotlib.pyplot as plt
from scipy.signal import gausspulse

def gausSin(F0, BW, Fs, phi=0):
    t = np.linspace(0,1,num=Fs,endpoint=False)
    var = 0.441 / BW
    a = 1 / (var*(2*np.pi)**0.5)
    y = np.sin(2*np.pi*F0*t+phi)*a*np.exp(-(t-0.5)**2/(2*var**2))
    # Normalize
    y /= np.absolute(y).max()
    return y

def next_power2(n):
    return 2.0**(2*n-1).bit_length()
    
## Generate dispersed signal

class dispSig():
    def __init__(self,fc,bw,dm,Fs,SNR,shiftlimit=200,window=1024,unit='hz'):
        """ fc: center frequency / bw: bandwidth / dm: Dispersion (pc*cm-3)/ Fs: Sampling Frequency\n
        SNR: Signal to Noise Ratio / timeL: Total time (s) / shiftlimit: remove shifts greater than limit\n
        window: length of window (?)"""
        
        self.unit=unit.lower()
        if unit in ['hz','khz','mhz','ghz']:
            if   unit ==  'hz': self.mag = 1.0e9
            elif unit == 'khz': self.mag = 1.0e6
            elif unit == 'mhz': self.mag = 1.0e3
            elif unit == 'ghz': self.mag = 1.0
        else:
            print 'Selected default unit: Hz'
            self.mag = 1.0
        self.fc = fc
        self.bw = bw
        self.dm = dm
        self.Fs = Fs
        self.SNR = SNR
        self.shiftlimit = shiftlimit
        self.window = window
    
    def makeSignal(self,n=0):
        if n==0:
            n = self.Fs
        bwf = 1.*self.bw/self.fc
        t = np.linspace(-.5,.5,num=n)
        self.signal = gausspulse(t,fc=self.fc,bw=bwf,bwr=-3)
    
    def createData(self):
        print 'Creating data...'
        noise = np.random.random(len(self.signal))
        self.data = self.SNR*self.signal+noise
        self.fpower = np.fft.rfft(self.data)
#        self.fpower = self.fpower[1:]
        print 'Data created.'
    
        # Dispersion
    def disperse(self):
        print 'Calculating dispersion...'
        f = np.fft.rfftfreq(int(self.Fs),d=1/self.Fs)
        pind = np.where(f > 0)
        self.fpower = self.fpower[pind]
        self.fPos = f[pind]
        self.sftT = 4.149e-3*self.dm*(self.mag/self.fPos)**2
        
        indic = np.where(self.sftT>self.shiftlimit)
        for i in indic:
            self.fpower[i] = 0.0
        phaseShft = np.exp(-1j*2*np.pi*self.fPos*self.sftT)
        fpowerShft = self.fpower*phaseShft
        self.TimeDisp = np.fft.irfft(fpowerShft)
        print 'Dispersion calculated.'
    
    # Create Base Bank data
    def createSignal(self):
    #    BBfreq = np.fft.fftfreq(window,d=1/Fs)
    ######################
    #### CREATE DATA
    ######################
        print 'Creating Base Bank data'
        self.BsBankD = np.abs(np.fft.rfft(self.TimeDisp[0:self.window]))
        BsBankL = len(self.TimeDisp)/self.window-1
        for i in xrange(BsBankL):
            self.BsBankD = np.vstack((np.abs(np.fft.rfft(
                           self.TimeDisp[i*self.window:(i+1)*self.window])),
                           self.BsBankD))
            # write percent complete        
            complete = 100.0*i/BsBankL
            sys.stdout.write("\r{0:.2f}".format(complete) + '%')
            sys.stdout.flush()
        print '\nData created'
        
    def showimg(self):
        plt.figure(1)
        plt.imshow(self.BsBankD.T,origin = "lower",interpolation='nearest',
                   extent=[0,self.timeL,0,self.fPos.max()],aspect='auto',cmap='gnuplot2')
        plt.ylim(self.fPos.min(),self.fPos.max())
        plt.ylabel("Frequency ({})".format(self.unit))
        plt.xlabel("Time (s)")
        #plt.figure(2)
        #plt.subplot(2,1,1)
        #plt.plot(data)
        #plt.subplot(2,1,2)
        #plt.plot(DispsTimeSris)
        #plt.figure(3)
        #plt.plot(sftT)
        plt.show()
    
    def auto(self,signal=True):
        if signal: self.makeSignal()
        self.createData()
        self.disperse()
        self.createSignal()
        self.showimg()
