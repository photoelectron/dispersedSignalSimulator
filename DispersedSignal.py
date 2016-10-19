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
    def __init__(self,fc,bw,dm,fs,SNR,shiftlimit=0,window=1024,unit='hz',timeL=1,zpad=False):
        """ fc: center frequency / bw: bandwidth / dm: Dispersion (pc*cm-3)/ fs: Sampling Frequency\n
        SNR: Signal to Noise Ratio / timeL: Total time (s) / \n
        shiftlimit: remove shifts greater than limit; if 0, then no removal\n window: length of window\n
        zpad: apply zero padding"""
        
        self.unit=unit.lower()
        if unit in ['hz','khz','mhz','ghz']:
            if   unit ==  'hz': self.mag = 1.0
            elif unit == 'khz': self.mag = 1.0e3
            elif unit == 'mhz': self.mag = 1.0e6
            elif unit == 'ghz': self.mag = 1.0e9
        else:
            print 'Selected default unit: Hz'
            self.mag = 1.0
        self.timeL = timeL # Currently not used
        self.fc = fc
        self.bw = bw
        self.dm = dm
        self.fs = fs
        self.snr = SNR
        self.shiftlimit = shiftlimit
        self.window = window
        self.zpad = zpad
    
    def makeSignal(self,n=0):
        """ Generates a Gaussian pulse given the object's parameters;
            n: number of samples; if 0, then sampling frequency is used."""
        if n==0:
            n = self.fs
        bwf = 1.*self.bw/self.fc
        t = np.linspace(-(self.timeL/2.)/self.mag,(self.timeL/2.)/self.mag,num=n)
        signal = gausspulse(t,fc=self.fc*self.mag,bw=bwf,bwr=-3)
        ### Zero padding
        if self.zpad:
            signal = np.pad(signal,int(self.fs*2),'constant')
        return signal
    
    def createData(self,signal):
        """ Adds noise to the signal and calculates FFT."""
        print 'Creating data...'
        noise = np.random.random(len(signal))
        tData = self.snr*signal+noise
        print 'Data created.'
        return tData
    
    def calcFFT(self,tData):
        """ Transforms time domain signal to freq domain"""
        fData = np.fft.rfft(tData)
        f = np.fft.rfftfreq(int(self.fs),d=1/self.fs)
        pind = np.where(f > 0)
        fData = fData[pind]
        self.fPos = f[pind]
        return fData
        
    # Dispersion
    def disperse(self,fData):
        """ Calculates time dispersion for each frequency, applies it
            in frequency domain, then calculates inverse FFT to obtain
            the dispersed signal in time."""
        print 'Calculating dispersion...'
        
        sftT = (4.149e-3/self.mag)*self.dm*(self.mag/self.fPos)**2
        if self.shiftlimit > 0:
            indic = np.where(sftT>self.shiftlimit)
            for i in indic:
                fData[i] = 0.0
        phaseShft = np.exp(-1j*2*np.pi*self.fPos*sftT)
        fpowerShft = fData*phaseShft
        dispData = np.fft.irfft(fpowerShft)
        print 'Dispersion calculated.'
        return dispData
    
    # Create Base Bank data
    def createSignal(self,dispData):
        print 'Creating Base Bank data'
        fbData = np.abs(np.fft.rfft(dispData[0:self.window]))
        fbDataL = len(dispData)/self.window-1
        for i in xrange(fbDataL):
            fbData = np.vstack((np.abs(np.fft.rfft(
                           dispData[i*self.window:(i+1)*self.window])),
                           fbData))
            # write percent complete        
            complete = 100.0*i/fbDataL
            sys.stdout.write("\r{0:.2f}".format(complete) + '%')
            sys.stdout.flush()
        print '\nData created'
        return fbData
        
    def showimg(self,img):
        """ Show the waterfall plot."""
        plt.figure(1)
        plt.imshow(img.T,origin = "lower",interpolation='nearest',
                   extent=[-self.timeL/2.,self.timeL/2.,0,self.fPos.max()],aspect='auto',cmap='gnuplot2')
        plt.ylim(self.fPos.min(),self.fPos.max())
        plt.ylabel("Frequency ({})".format(self.unit))
        plt.xlabel("Time (s)")
        plt.show()
    
#    def auto(self,signal=True):
#        """ Runs all functions to create and show a dispersed signal.
#            If signal is False, it does not generate a new Gaussian pulse."""
#        if signal: self.makeSignal()
#        self.createData()
#        self.disperse()
#        self.createSignal()
#        self.showimg()
    
    def auto2(self,signal=True):
        """ Runs all functions to create and show a dispersed signal.
            If signal is False, it does not generate a new Gaussian pulse."""
        if signal: 
            signal = self.makeSignal()
        tData = self.createData(signal)
        fData = self.calcFFT(tData)
        dispData = self.disperse(fData)
        fbData = self.createSignal(dispData)
        self.showimg(fbData)
