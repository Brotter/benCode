import scipy
import numpy as np
import pylab as lab

import tfUtils as tf


def normalize(inGraphY):

    inGraphY /= np.std(inGraphY)
    inGraphY /= np.sqrt(len(inGraphY))

    return inGraphY


def makeNoise(n, f_lo=0.15, f_hi=1.2, plot=False):

    """
    Kieth's noise generator

    Sampling rate in GHz
    Low and high frequency boundaries in GHz
    """

    t_total = n / 10 # time of trace in ns
    if n % 2 == 0:
        index_nyquist = n / 2
    else:
        index_nyquist = (n - 1) / 2
    f = np.array(range(1, index_nyquist + 1),dtype=float) / t_total # Frequency in GHz
    cut = np.logical_and(f >= f_lo, f <= f_hi)
    phi = np.random.uniform(0, 2. * np.pi, index_nyquist)
    amplitude = scipy.stats.rayleigh.rvs(size=index_nyquist)
    a_postive_frequency = cut * amplitude * (np.cos(phi) + np.complex(0, 1) * np.sin(phi))
    #a = np.zeros(n) + np.zeros(n) * np.complex(0, 1)
    a = np.ones(n) * np.complex(0., 0.)
    a[1: index_nyquist + 1] = a_postive_frequency
    
    y = np.real(np.fft.ifft(a))
    y /= np.std(y)
    y /= np.sqrt(n) #need this to make the autocorrelation equal to one

    return y




def compNoise(n,numPts = 300):

    """

    Generate an expectation of the peak correlation value from just noise

    """


    peaks = []

    firstNoise = makeNoise(numPts)
    for i in range(0,n):
        corr = tf.correlation(firstNoise,makeNoise(numPts))
        max = np.max(corr)
        min = np.min(corr)
        peaks.append(np.max([max,-min]))

    lab.hist(peaks,np.arange(0,100)/100.)
    lab.show()


def addNoise(sigWaveX,sigWaveY,snr,showPlots=False):

    length = len(sigWaveX)

    sigWaveY /= np.std(sigWaveY)
    sigWaveY /= np.sqrt(length)
    sigWaveY *= snr

    noiseX = np.arange(0,length)*0.1
    noiseY = makeNoise(length) #normalized to one already

    noiseF,noiseFFT = tf.genFFT(noiseX,noiseY)
    sigF,sigFFT = tf.genFFT(sigWaveX,sigWaveY)

    convolved = sigFFT+noiseFFT

    outX,outY = tf.genTimeSeries(sigF,convolved)

    outY /= np.std(outY)
    outY /= np.sqrt(length)


    if showPlots:
        fig,ax = lab.subplots(3,2)
        ax[0][0].plot(sigWaveX,sigWaveY,color="blue",label="signal")
        ax[0][0].legend()
        ax[2][0].plot(noiseX,noiseY,color="red",label="noise")
        ax[2][0].legend()
        ax[1][0].plot(outX,outY,color="purple",label="convolution")
        ax[1][0].legend()

        ax[0][1].plot(sigF,tf.calcLogMag(sigF,sigFFT),color="blue",label="signal")
        ax[0][1].legend()
        ax[2][1].plot(noiseF,tf.calcLogMag(noiseF,noiseFFT),color="red",label="noise")
        ax[2][1].set_ylim([-30,30])
        ax[2][1].legend()
        ax[1][1].plot(sigF,tf.calcLogMag(sigF,convolved),color="purple",label="convolution")
        ax[1][1].set_ylim([-30,30])
        ax[1][1].legend()
        fig.show()
        

    return outX,outY
    

def genHistAtSnr(n=100,snr=1,showPlots=False):

    templates = np.loadtxt("templates.txt").T

    corrs = []

    for trial in range(0,n):
        convX,convY = addNoise(templates[0],templates[1],snr)
        corr = np.max(tf.correlation(normalize(convY),normalize(templates[1])))
        corrs.append(corr)

    if showPlots:
        fig,ax = lab.subplots()
        ax.hist(corrs,np.arange(-0.1,1.1,0.001))
        fig.show()
    

    return corrs
        

def valueVsSnr(n=100):

    snrs = []
    means = []
    sigmas = []

    for snr in np.arange(0,5,0.01)+0.01:
        hist = genHistAtSnr(n=n,snr=snr)
        mean = np.mean(hist)
        sigma = np.std(hist)
        print snr,mean,sigma
        snrs.append(snr)
        means.append(mean)
        sigmas.append(sigma)
        


    fig,ax = lab.subplots(2)
    ax[0].set_title("Template value vs SNR")
    ax[0].plot(snrs,means)
    ax[0].set_ylabel("Mean Template Value")
    ax[1].plot(snrs,sigmas)
    ax[1].set_ylabel("Sigma")
    ax[1].set_xlabel("SNR")
    fig.show()


    return snrs,means,sigmas


def genHistAtWindow(nPts,nTrials=1000,snr=2,showPlots=False,templateNum=1,offset=120):
    templates = np.loadtxt("templates.txt").T

    initLen = len(templates[0])-offset
    changeInPts = nPts-initLen

    if (changeInPts > 0):
        tempWindX,tempWindY = tf.zeroPadEnd(templates[0],templates[1],changeInPts)
    elif (changeInPts < 0):
        tempWindX = templates[0][offset:nPts+offset]
        tempWindY = templates[templateNum][offset:nPts+offset]
    else:
        tempWindX = templates[0][offset:]
        tempWindY = templates[templateNum][offset:]

    print len(tempWindX)

    corrs = []

    if showPlots:
        fig0,ax0 = lab.subplots()
        ax0.plot(tempWindX,tempWindY)
        fig0.show()

    for trial in range(0,nTrials):
        convX,convY = addNoise(tempWindX,tempWindY,snr,showPlots=(trial==0 and showPlots))
        corr = np.max(tf.correlation(normalize(convY),normalize(tempWindY)))
        corrs.append(corr)

    if showPlots:
        fig,ax = lab.subplots()
        ax.hist(corrs,np.arange(-0.1,1.1,0.001))
        fig.show()
    

    return corrs
        



def valueVsWindowLen(snr=2,templateNum=1,offset=0):

    """
    Also important is to see if I should be windowing these things
    """

    windows=[]
    means=[]
    sigmas=[]
    for window in np.arange(200,5000,50):
        hist = genHistAtWindow(window,snr=snr,templateNum=templateNum,offset=offset)
        windows.append(window)
        means.append(np.mean(hist))
        sigmas.append(np.std(hist))
        print window,means[-1],sigmas[-1]



    fig,ax = lab.subplots(2)
    ax[0].plot(windows,means)
    ax[1].plot(windows,sigmas)
    fig.show()

    return windows,means,sigmas
                     
