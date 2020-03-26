#QAM_Modulator

import matplotlib.pyplot as plt
from numpy.random import*
from numpy import*
from pylab import*

def rcosfilter(N, beta, Ts, Fs):
    t = (np.arange(N) - N / 2) / Fs
    return (np.where(np.abs(2*t) == Ts / beta,
        np.pi / 4 * np.sinc(t/Ts),
        np.sinc(t/Ts) * np.cos(np.pi*beta*t/Ts) / (1 - (2*beta*t/Ts) ** 2)),t)

#Symbol rate
Sr = 1e6

#Upsampling factor
K = 8

#Number of symbols
Ns = 128

#SNR
snr = 35

#Scattering coeficient of shaping filter
shapeBeta = 0.35

#Span of filter in symbols
shapeSpan =  8

#Number of bits per symbol
Bs = 4

#Constellation
QAM16 = [-1, -0.33, 0.33, 1]

#Intermediate frequency
IF = 2e6

#Binary data
data = randint (0,2,Ns*Bs)

#Data
dataI = data[0::2]
dataQ = data[1::2]

#Mapper and Upsampler
upI = [0]*K*Ns
upQ = [0]*K*Ns

for i in range(0,len(dataI),2):
	upI[(i/2)*K] = QAM16[2*dataI[i]+ dataI[i+1]]
	upQ[(i/2)*K] = QAM16[2*dataQ[i]+ dataQ[i+1]]

# Raised cosine filter 
h,t = rcosfilter (shapeSpan*K, shapeBeta, 1/Sr, Sr*K) 

shapedI = np.convolve(upI, h)
shapedQ = np.convolve(upQ, h)

shapedI = shapedI[shapeSpan*K/2:1+len(shapedI)-shapeSpan*K/2]
shapedQ = shapedQ[shapeSpan*K/2:1+len(shapedQ)-shapeSpan*K/2]

#NCO
t = linspace(0,Ns/Sr,Ns*K)
loCos = cos(2*pi*IF*t) 
loSin = sin(2*pi*IF*t)

#Mixers
mixI = np.multiply(loCos,shapedI)
mixQ = np.multiply(loSin,shapedQ)

#Combiner
signal = mixI + mixQ

#Plots
t = linspace(0,Ns/Sr,Ns*Bs)
figure('Data stream')
subplot(3,1,1)
stem(t,data)

t = linspace(0,Ns/Sr,Ns*Bs/2)
subplot(3,1,2)
stem(t,dataI)

subplot(3,1,3)
stem(t,dataQ)

figure('Saidas do Upsampler')
t = linspace(0,Ns/Sr,Ns*K)
subplot(2,1,1)
stem(t,upI)
subplot(2,1,2)
stem(t,upQ)

figure('Saidas do filtro')
subplot(3,1,1)
plot(t, shapedI)
subplot(3,1,2)
plot(t, shapedQ)
subplot(3,1,3)
plot(shapedI,shapedQ)

figure('Saida do Combiner')
plot(t,signal)


figure('Resposta ao impulso do filtro')
t = linspace(0,shapeSpan/Sr,shapeSpan*K) - shapeSpan/(2*Sr)
stem(t,h)



signal = upI
# Plotando a FFT
FFT = np.fft.fft(signal)
DC  = 20*log((2.0/len(signal)) * abs(FFT[0]))
amp = 20*log((2.0/len(signal)) * abs(FFT[1:len(signal)/2+1:1]))
amp = np.insert(amp,0,DC)
f = arange (0,len(signal)/2+1, 1) * Sr * K / (len(signal)*1000)

figure('Espectro')
plot(f,amp)



show()
