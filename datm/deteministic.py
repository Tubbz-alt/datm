#!/usr/bin/python

import numpy as np;
from numpy.fft import fft as FFT;
from numpy.fft import ifft as IFFT;
from numpy.fft import fftfreq as FREQS;
from os import listdir;
from os.path import isfile, join;
import fnmatch;

def gauss(x,x0,w):
	return np.exp(-1*np.power((x-x0)/w,int(2)));

dirname = '/data/raw/';
procdir = '/datm/data/processed/';
filenames = [f for f in listdir(dirname) if isfile(join(dirname, f)) ];
print(len(filenames));
filebase = 'noetalon_1200_interference.out.'
#filebase = 'nfibers109_617_interference.out.'
#/data/raw/nfibers109_617_interference.out.9990
for filename in filenames:
	if fnmatch.fnmatch(filename, filebase + '*'):
		print(filename);
		data = np.loadtxt(dirname+filename,dtype=float);
		break;
print('nfibers: data.shape[0] = ',data.shape[0]);
print('npixels: data.shape[1] = ',data.shape[1]);

fiberids = np.arange(data.shape[0],dtype=float);
print(fiberids);

freqs = np.zeros(data.shape,dtype=float);
nfiles=1000;
catindsout=np.zeros((0,6),dtype=float);
indsout = np.zeros((data.shape[0],6),dtype=float);
ntally=int(0);
overwidth=.02/2.;# .01 for noetalon
atwidth=.1/2.; # .1 for noetalon
for filename in filenames:
	if ntally > nfiles:
		break;
	if fnmatch.fnmatch(filename, filebase+'*'):
		#print(filename);
		ntally +=1;
		f = open(dirname+filename, 'r');
		s = f.readline();
		string,val = [splits for splits in s.split("\t") if splits is not ""];
		delay = float(val);
		data = np.loadtxt(dirname+filename,dtype=float);
		dataFT=FFT(data,axis=1);
		#dataFTabs=np.abs(dataFT);
		dataFTfilt = np.zeros(dataFT.shape,dtype=complex);
		dataFToverfilt = np.zeros(dataFT.shape,dtype=complex);
		f = list(FREQS (dataFT.shape[1])) * dataFT.shape[0];
		freqs = np.reshape(f,dataFT.shape);
		dataFT *= 1j*freqs;
		dataFToverfilt = dataFT*gauss(freqs,0,overwidth);
		dataFTfilt = dataFT*gauss(freqs,0,atwidth);

		dataBack = IFFT(dataFTfilt,axis=1);
		dataOverBack = IFFT(dataFToverfilt,axis=1);
		out = np.real(dataBack)*np.abs(dataOverBack);
		maxs = np.argmax(out,axis=1);
		mins = np.argmin(out,axis=1);
		maxvals = np.max(out,axis=1);
		minvals = np.min(out,axis=1);
		indsout[:,0] = int(delay) ; # in femtoseconds
		indsout[:,1] = fiberids; 
		indsout[:,2] = maxs;
		indsout[:,3] = mins;
		indsout[:,4] = maxvals;
		indsout[:,5] = minvals;
		catindsout = np.row_stack((catindsout,indsout));


		#np.savetxt(procdir+filename+'.fftabs',dataFTabs,fmt='%.3e');
		if (ntally%50==0):
			np.savetxt(procdir+filename+'.back',out,fmt='%.3e');
			np.savetxt(procdir+filebase+'allinds.out',catindsout,fmt='%.3f');

np.savetxt(procdir+filebase+'allinds.out',catindsout,fmt='%.3f');

