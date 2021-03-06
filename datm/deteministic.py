#!/usr/bin/python

import numpy as np;
from numpy.fft import fft as FFT;
from numpy.fft import ifft as IFFT;
from numpy.fft import fftfreq as FREQS;
from os import listdir;
from os.path import isfile, join;
import fnmatch;
import pdb as debug;

def gauss(x,x0,w):
	return np.exp(-1*np.power((x-x0)/w,int(2)));

def courseTime(x):
	sz = x.shape[1];
	filter = list( gauss(np.arange(sz),sz//2,sz//8)) * x.shape[0];
	filter = np.roll(sz//2,axis=1);
	filterFT = FFT(filter,axis=1);
	xFT = FFT(x,axis=1);
	yFT = xFT*filterFT;
	y = IFFT(yFT,axis=1);
	i = argmax(y,axis=1);
	return (i,y);


dirname = '/data/raw/';
procdir = '/datm/data/processed/';
filenames = [f for f in listdir(dirname) if isfile(join(dirname, f)) ];
print(len(filenames));
filebase = 'uint12_617_interference.out.'
#filebase = 'noetalon_1200_interference.out.'
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

ncols=int(10);

freqs = np.zeros(data.shape,dtype=float);
nfiles=1000;
catindsout=np.zeros((0,7),dtype=float);
indsout = np.zeros((data.shape[0],7),dtype=float);
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
		listofstrings = [strings for strings in filename.split(".") if strings is not ""];
		imagenum = listofstrings[-1];
		delay = float(val);
		data = np.loadtxt(dirname+filename,dtype=float);
		coarseinds,coarseconv = courseTime(np.copy(data));
		print(coarseinds[:10]);
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
		#maxs = np.argmax(out,axis=1);
		#mins = np.argmin(out,axis=1);
		sortinds=np.argsort(out,axis=1);
		#print(out.shape[0]);
		maxs_f = np.zeros(out.shape[0],dtype=float);
		mins_f = np.zeros(out.shape[0],dtype=float);
		maxvals = np.zeros(out.shape[0],dtype=float);
		minvals = np.zeros(out.shape[0],dtype=float);
		for i in range(out.shape[0]):
			maxvals[i] = np.sum(out[i,sortinds[i,-ncols:]]);
			minvals[i] = np.sum(out[i,sortinds[i,:ncols]]);
			if (np.abs(maxvals[i])+np.abs(minvals[i])) < 1:
				maxvals[i] = 0;
				minvals[i] = 0;
				continue;
			#print(out[i,sortinds[i,:3]],out[i,sortinds[i,-3:]]);
			maxs_f[i] = np.average(sortinds[i,-ncols:],weights=np.abs(out[i,sortinds[i,-ncols:]]));
			mins_f[i] = np.average(sortinds[i,:ncols],weights=np.abs(out[i,sortinds[i,:ncols]]));
			#print(maxs_f[i],mins_f[i]);
#### AHHHH I FUCKING HATE PYTHON !!!!!!!!!!!!
		#sortinds=np.unravel_index(np.argsort(out,axis=1),out.shape);
		#print(out[sortinds[:,-ncols:]]);
		#for i in range(4)
		#break;
		#maxvals = np.max(out,axis=1);
		#minvals = np.min(out,axis=1);
		indsout[:,0] = int(delay) ; # in femtoseconds
		#debug.set_trace();
		indsout[:,1] = int(imagenum) ; 
		indsout[:,2] = fiberids; 
		indsout[:,3] = maxs_f;
		indsout[:,4] = mins_f;
		indsout[:,5] = maxvals;
		indsout[:,6] = minvals;
		catindsout = np.row_stack((catindsout,indsout));


		#np.savetxt(procdir+filename+'.fftabs',dataFTabs,fmt='%.3e');
		if (ntally%50==0):
			np.savetxt(procdir+filename+'.back',out,fmt='%.3e');
			np.savetxt(procdir+filebase+'allinds.out',catindsout,fmt='%.3f');

np.savetxt(procdir+filebase+'allinds.out',catindsout,fmt='%.3f');

