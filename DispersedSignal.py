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
    def __init__(self,fc,bw,dm,fs,SNR,shiftlimit=200,window=1024,unit='hz'):
        """ fc: center frequency / bw: bandwidth / dm: Dispersion (pc*cm-3)/ fs: Sampling Frequency\n
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
        self.timeL = 1
        self.fc = fc
        self.bw = bw
        self.dm = dm
        self.fs = fs
        self.snr = SNR
        self.shiftlimit = shiftlimit
        self.window = window
    
    def makeSignal(self,n=0):
        """ Generates a Gaussian pulse given the object's parameters;
            n: number of samples; if 0, then sampling frequency is used."""
        if n==0:
            n = self.fs
        bwf = 1.*self.bw/self.fc
        t = np.linspace(-.1,.1,num=n)
        self.signal = gausspulse(t,fc=self.fc,bw=bwf,bwr=-3)
        ### Zero padding
        self.signal = np.pad(self.signal,int(self.fs*2),'constant')
    
    def createData(self):
        """ Adds noise to the signal and calculates FFT."""
        print 'Creating data...'
        noise = np.random.random(len(self.signal))
        self.data = self.snr*self.signal+noise
        self.fpower = np.fft.rfft(self.data)
#        self.fpower = self.fpower[1:]
        print 'Data created.'
    
    # Dispersion
    def disperse(self):
        """ Calculates time dispersion for each frequency, applies it
            in frequency domain, then calculates inverse FFT to obtain
            the dispersed signal in time."""
        print 'Calculating dispersion...'
        f = np.fft.rfftfreq(int(self.fs),d=1/self.fs)
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
        """ Show the waterfall plot."""
        plt.figure(1)
        plt.imshow(self.BsBankD.T,origin = "lower",interpolation='nearest',
                   extent=[0,self.timeL,0,self.fPos.max()],aspect='auto',cmap='gnuplot2')
        plt.ylim(self.fPos.min(),self.fPos.max())
        plt.ylabel("Frequency ({})".format(self.unit))
        plt.xlabel("Time (s)")
        plt.show()
    
    def auto(self,signal=True):
        """ Runs all functions to create and show a dispersed signal.
            If signal is False, it does not generate a new Gaussian pulse."""
        if signal: self.makeSignal()
        self.createData()
        self.disperse()
        self.createSignal()
        self.showimg()
