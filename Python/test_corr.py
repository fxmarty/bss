Fs <- 8000
dt <- 1/Fs
f1 <- 100
tdelay <- 0.625e-03
t3 <- seq(0,1-dt,dt)
x3 <- cos(2*pi*f1*t3)
x4 <- cos(2*pi*f1*(t3-tdelay))

xcorr_31956 <- function(x,y,normalize = FALSE)
{
  xfft <- fft(x, 4*length(x))
  yfft <- fft(y, 4*length(x))
  
  R <- xfft*Conj(yfft);
  if (normalize)
  {
    R <- R/abs(R)
  }
  c <- fft(R, inverse=TRUE);
  
  return(c)
}

xc <- xcorr_31956(x3,x4, FALSE)
xc_phat <- xcorr_31956(x3,x4, TRUE)

par(mfrow=c(2,1))
plot(seq(0,length(xc)-1),Re(xc), type="l", xlim=c(0,20), col="blue", lwd=2)
ix <- which.max(abs(xc))
points(ix-1,abs(xc[ix]), col="red", lwd=5); 
title('Standard CCF')

plot(seq(0,length(xc)-1),Re(xc_phat), type="l", xlim=c(0,20), col="blue", lwd=2)
ix_phat <- which.max(abs(xc_phat))
points(ix_phat-1,abs(xc_phat[ix_phat]), col="red", lwd=5); 
title('Generalized')

print(paste("Delay is",tdelay*Fs), quote = FALSE)
print(paste("Estimate is",which.max(abs(xc)) - 1), quote = FALSE)

##
    xfft = np.fft.fft(x3, 4*len(x3))
    yfft = np.fft.fft(x4, 4*len(x3))
    R = xfft * np.conjugate(yfft)
    if normalize:
        R = R/np.abs(R)
    c = np.fft.ifft(R)
    return c
##
plt.plot(t3,x3)
plt.grid()
plt.show()
    
##
Fs = 8000
dt = 1/Fs
f1 = 100
tdelay = 20/8000
t3 = np.arange(0,1,dt)
x3 = np.cos(2*np.pi*f1*t3) #+ np.random.normal(0, 0.05, size=len(t3))

x4 = np.cos(2*np.pi*f1*(t3-tdelay)) #+ np.random.normal(0, 0.05, size=len(t3))

#x3 = np.cos(2*np.pi*f1*t3)
#x4 = np.cos(2*np.pi*f1*(t3-tdelay))

def xcorr_31956(x,y,normalize):
    #xfft = np.fft.fft(x, 4*len(x))
    #yfft = np.fft.fft(y, 4*len(x))
    xfft = np.fft.fft(x, 4*len(x))
    yfft = np.fft.fft(y, 4*len(x))
    R = xfft * np.conjugate(yfft)
    if normalize:
        R = R/(np.abs(R))
    c = np.fft.ifft(R)
    return c

xc = xcorr_31956(x3,x4,False)
xc_phat = xcorr_31956(x3,x4,True)

abs = np.arange(-len(x3),len(x3),1)


#par(mfrow=c(2,1))
#plot(seq(0,length(xc)-1),Re(xc), type="l", xlim=c(0,20), col="blue", lwd=2)

bon1 = np.fft.fftshift(np.abs(xc))[len(xc)//2 - len(x3):len(xc)//2 + len(x3)]
bon2 = np.fft.fftshift(np.abs(xc_phat))[len(xc)//2 - len(x3):len(xc)//2 + len(x3)]
bon2[len(x3)] = 0

plt.subplot(2, 1, 1)
plt.plot(abs, bon1)
plt.title('A tale of 2 subplots')
plt.xlim(xmax=100)  
plt.xlim(xmin=-100)
plt.grid()


plt.subplot(2, 1, 2)
#plt.plot(abs, np.fft.fftshift(np.real(xc_phat)))
plt.plot(abs, bon2)
plt.xlim(xmax=100)  
plt.xlim(xmin=-100)

plt.grid()
plt.show()


##
xfft = np.fft.fft(x3, 4*len(x3))
yfft = np.fft.fft(x4, 4*len(x3))
R = xfft * np.conjugate(yfft)
R = R/np.abs(R)
c = np.fft.ifft(R)
c = np.fft.fftshift(c)

##
for i in range(0,-30,-1):
    print(i,bon2[8000+i])

##

plt.plot(np.linspace(0,1,len(res1)),res1_red)
plt.grid()
plt.show()
##
def xcorr_31956(x,y,normalize):
    #xfft = np.fft.fft(x, 4*len(x))
    #yfft = np.fft.fft(y, 4*len(x))
    xfft = np.fft.fft(x, 4*len(x))
    yfft = np.fft.fft(y, 4*len(x))
    R = xfft * np.conjugate(yfft)
    if normalize:
        R = R/np.abs(R)
    c = np.fft.ifft(R)
    return c

#res1_red = m1[0:64000]
#res2_red = m2[0:64000]

bruit = np.random.normal(0, 0.005, size=len(res1)+40)
#+ 0.4*bruit[20:len(res1)+20]
#+ bruit[0:len(res1)]

#res1_red = res1 + 0.4*bruit[len(res1)] + np.random.normal(0, 0.05, size=len(res1))
#res2_red = 0.4*res2 + bruit[24:len(res1)+24] + np.random.normal(0, 0.05, size=len(res1))

res1_red = res1
res2_red = res2


xc = xcorr_31956(res1_red,res2_red,False)
xc_phat = xcorr_31956(res1_red,res2_red,True)

abs = np.arange(-len(res1_red),len(res1_red),1)


bon1 = np.fft.fftshift(np.real(xc))[len(xc)//2 - len(res1_red):len(xc)//2 + len(res1_red)]
bon2 = np.fft.fftshift(np.real(xc_phat))[len(xc)//2 - len(res1_red):len(xc)//2 + len(res1_red)]
bon2[len(res1_red)] = 0

plt.subplot(2, 1, 1)
plt.plot(abs, bon1)
plt.title('A tale of 2 subplots')
plt.xlim(xmax=100)  
plt.xlim(xmin=-100)
plt.grid()


plt.subplot(2, 1, 2)
#plt.plot(abs, np.fft.fftshift(np.real(xc_phat)))
plt.plot(abs, bon2)
plt.xlim(xmax=100)  
plt.xlim(xmin=-100)

plt.grid()
plt.show()
