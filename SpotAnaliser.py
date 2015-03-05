import matplotlib as mpl
import matplotlib.pyplot as plt
import prettyplotlib as ppl
import pims
import numpy as np
import scipy as sp
from scipy import ndimage
from scipy.interpolate import UnivariateSpline
from scipy.stats import linregress
mpl.rc('figure',  figsize=(20, 20))
mpl.rc('image', cmap='gray', interpolation='none')

import os
import gc

class Spot_Analysis:
	def __init__(self, filename):
		#################################
		#Loading the data and starting the preprocessing
		#################################
		self.filename = filename
		self.img = pims.TiffStack(filename)
		self.time_pool, self.binned = self.Pool_Time(self.img)
		#################################

		#################################
		#Definition of constants
		#################################
		self.deltax = 0.1 						#Defines the delta x step for derivation etc
		self.mid = len(self.time_pool)/2 			#Define Middle pixel column
		self.pix_col = range(self.mid-1,self.mid+2)		#Defines the X pixel to use for speed analysis
		self.timeinterval = 2 					# Time between Frames in second 
		self.speed_threshold = 40 				#Speed Threshold, if an object moves faster than this, it is excluded !
		self.timethreshold = 4 *self.timeinterval*(1/self.deltax) #Defined how far appart can peak be (in second)
		self.Time_Error = 6 / self.timeinterval #Allowed variance in second for a peak to be called true when found in vertical axis and checked VS time  
		self.diameter = 850 					#Diameter of the bacteria in nm
		self.pix_size = 65 						#Size in nm of the pixel (here CMOS camera with 100x objective = 65nm)
		#################################

		#################################
		#Definition of Place Holders
		#################################
		self.all_v = []
		self.All_hit = []
		self.All_time = {}
		self.Real_hit = []
		#################################

	def bin_img(self, img, Sum=True, width=1, height=1):
		"Bin a 2D picture with binning of size Width and height. If Sum==True, then sums to pixels, otherwise, average them"
		res = []
		if width>1:
			for im in img.T:
				data = im
				if Sum:
					res.append(data[:(data.size // width) * width].reshape(-1, width).sum(axis=1))
				else:
					res.append(data[:(data.size // width) * width].reshape(-1, width).mean(axis=1))
			res = np.array(res)
		else:
			res = img.T
		if height>1:
			res2 = []
			for im in res:
			   data = im
			   if Sum:
				   res2.append(data[:(data.size // width) * width].reshape(-1, width).sum(axis=1))
			   else:
				   res2.append(data[:(data.size // width) * width].reshape(-1, width).mean(axis=1))
			res2 = np.array(res2).T
			res = res2
		return res.T

	def bin_stack(self, stack, width=1, Sum=True):
		"Perform binning on a Tiff Stack (multipage)"
		binned_stack = []
		for i in range(len(stack)):
			binned_stack.append(self.bin_img(stack[i], width=width, Sum=Sum))
		return binned_stack

	def find_nearest(self, array, value):
		"Finds the nearest element in the list and return its indice"
		idx = (np.abs(array-value)).argmin()
		return array[idx]

	def Pool_Time(self, img, Bin=False):
		"Transforms the pictures into a time matrice"
		if Bin:
			binned = self.bin_stack(img)
		else:
			binned = img
		time_pool = []

		for i in range(len(binned[0].T)):
			tmp = []
			for im in binned:
				tmp.append(im.T[i])
			time_pool.append(tmp)
		time_pool = np.array(time_pool)
		time_pool = time_pool - time_pool.min()
		return time_pool, binned

	def Fit_Affine(self, Y):
		"Fit affine lines through the data to correct for potential Decay"
		X = np.linspace(0,len(Y)*self.timeinterval,len(Y))
		s = linregress(X,Y)
		return s[0],s[1]

	def Signal_Process(self, Y, Gaussian=True, Gaussian_sigma=1, affine_correct=False, Log=True, MinRemove=True):
		if Gaussian:
			Y = ndimage.filters.gaussian_filter1d(Y, sigma=Gaussian_sigma, mode='reflect')
		if Log:
			Y = np.log(Y)
		if affine_correct:
			tY = []
			a,b = Fit_Affine(Y) 
			for i in range(len(X)):
				if (a*X[i] + b) > 0:
					tY.append(Y[i]-(a*X[i] + b))
				else:
					tY.append(Y[i])
			Y = np.array(tY)
		if MinRemove:
			Y = Y-Y.min()
		return Y

	def Peak_detect(self, Y, line=0, y_threshold=0.4, dy_threshold = 0.001):
		#Process Signal
		Y = self.Signal_Process(Y)

		#Create X vectors
		X = np.linspace(0,len(Y)*self.timeinterval,len(Y))
		xs = np.arange(0, len(Y)*self.timeinterval, self.deltax)

		#Fit Spline model to the data
		spl = UnivariateSpline(X, Y, k=5, s=0)

		#Take the Spline derivative to look for d(y) = 0 and dd(y) <= 0
		dspl = spl.derivative()   		#first derivative
		ddspl = dspl.derivative() 		#second derivative

		#generate Smooth data from the curve spline fitting 
		y = spl(xs)						
		dy = dspl(xs)
		ddy = ddspl(xs)

		#Find the roots, dy == 0 and ddy <= 0
		roots = []
		for i in range(len(dy)):
			if -dy_threshold < dy[i] and dy[i] < dy_threshold and y[i] > y_threshold and ddy[i] <= 0:
				roots.append(i)

		#Filter out the points that are next to each other  
		roots2 = []
		for i in range(len(roots)):
			if i < len(roots)-1:
				if not (roots[i+1] - roots[i]) < 60:
					roots2.append(roots[i])
			else:
				roots2.append(roots[i])
		roots = roots2
		return roots, X, xs, Y, y, dspl(xs), ddspl(xs)

	def Filter_Neighbor_Hit(self,roots_all, xs):
		hitlist = []
		roots_ok = []
		for r in roots_all[1]:
				for r0 in roots_all[0]:
					for r2 in roots_all[2]:
						if abs(r - r2) < self.timethreshold and abs(r - r0) < self.timethreshold and abs(r2 - r0) < 2*self.timethreshold:  
							if r0 < r < r2 or r2 < r < r0:	
								roots_ok.append((r0,'b'))
								roots_ok.append((r,'g'))
								roots_ok.append((r2,'r'))
							else:
								roots_ok.append((r0,'k'))
								roots_ok.append((r,'k'))
								roots_ok.append((r2,'k'))
							t = np.mean((abs(xs[r]-xs[r0]),abs(xs[r]-xs[r2]),abs(xs[r0]-xs[r2])/2))
							v = 65.0/t
							if v < self.speed_threshold:
								self.all_v.append(v)
							hitlist.append(int(round(xs[r],0)))
		return roots_ok, hitlist


	def Main(self, makeplot=False):
		time_hit = []
		self.bact_lenght = len(self.time_pool[0].T)
		self.bact_lenght_nm = self.bact_lenght*self.pix_size
		print "Bacteria Lenght:", len(self.time_pool[0].T)
		for line in range(len(self.time_pool[0].T)):
			# print "line is :",line
			roots_all = []
			for tp in self.pix_col:
				Y = self.time_pool[tp].T[line]
				roots, X, xs, Y, y, dy, ddy = self.Peak_detect(Y)
				roots_all.append(roots)
			roots_ok, hitlist= self.Filter_Neighbor_Hit(roots_all, xs)
			self.All_hit.append(np.unique(hitlist)/2)
			#Plot vertical line on the dectected peaks
			if makeplot:
				fig = plt.figure()
				plt.plot(X, Y ,label='pix%s'%tp)
				plt.plot(xs, y,label='pix%s'%tp)
				plt.plot(xs, dy, 'orange')
				plt.plot(xs, ddy, 'c')
				plt.plot(xs, [0]*len(xs), 'y')
			for i in roots_ok:
				if makeplot:
					plt.plot((xs[i[0]],xs[i[0]]), (0,3), i[1])
				if i[1] == 'g':
					time_hit.append(self.find_nearest(X,xs[i[0]]))
			if makeplot:
				plt.legend()
				fig.savefig(self.filename.split('.')[0] + '_time_peak_detect' + str(line) + '.png')
		time_hit = np.unique(time_hit)
		time_hit = [np.where(X == i)[0][0] for i in time_hit]
		self.speed = np.mean(self.all_v)
		print "mean Speed", self.speed
		for time in time_hit:
			# print "time is :",time
			Y = self.Signal_Process(self.time_pool[self.mid][time], Gaussian=False)
			roots, X, xs, Y, y, dy, ddy = self.Peak_detect(Y)			
			self.All_time[time] = [int(round(self.find_nearest(X,xs[i]),0)) for i in roots]
			#Plot vertical line on the dectected peaks
			if makeplot:
				fig = plt.figure()
				plt.plot(X, Y)
				plt.plot(xs,y)
				plt.plot(xs,dy,'orange')
				plt.plot(xs,ddy,'c')
				plt.plot(xs, [0]*len(xs), 'y')
				for i in roots_ok:
					plt.plot((xs[i],xs[i]), (0,3))
				plt.legend()
				fig.savefig(self.filename.split('.')[0] + '_vertical_peak_detect' + str(line) + '.png')
		

		i = 0
		for t in self.All_time.keys():
			for line in self.All_time[t]:
				test = False
				line = line-1
				for l in range(line-1,line+2):
					try:
						n = self.find_nearest(self.All_hit[line], t)
						# print "ERROR: ",line,n,t
						i+=1
						if abs(n-t) < self.Time_Error:
							test = True
					except:
						pass
				if test:
					self.Real_hit.append((line,t))

		Tn = {}
		lines, times = np.array(self.Real_hit)[:,0],np.array(self.Real_hit)[:,1]
		for i in range(len(times)):
			if Tn.has_key(times[i]):
				Tn[times[i]].append(lines[i])
			else:
				Tn[times[i]] = [lines[i]]
		T = np.unique(times)
		res = np.zeros(T.max())

		for t in range(len(T)):
			c = 0
			for l in Tn[T[t]]:
				if t < len(T)-2:
					if l not in Tn[T[t+2]] and l not in Tn[T[t+1]]:
						c+=1
				else:
					c+=1
			res[t] = c

		self.revolution_window = int((self.diameter*np.pi)/self.speed)
		windowed = []
		for i in range(len(res)-self.revolution_window):
			windowed.append(np.sum(res[i:i+self.revolution_window]))
		self.meanMREB, self.stdMREB = np.mean(windowed), np.std(windowed)
		print "Number of Events per revolution window:"
		print "Mean:",self.meanMREB, "Std:",self.stdMREB
		
		if makeplot:
			fig = plt.figure()
			plt.plot(range(len(windowed)),windowed)
			plt.title('Number of MReB throught the sliding window, Mean=%s, Std=%s'%(self.meanMREB,self.stdMREB))
			fig.savefig(self.filename.split('.')[0] + 'MReB_Detection.png')
			fig.clf()
			for i in range(len(self.binned)):
				fig = plt.figure()
				plt.imshow(self.binned[i])
				tt = np.where(times == i)[0]
				plt.plot(0,0,'k.')
				for t in tt:
					plt.plot(mid,lines[t],'ro')
				c = str(i)
				while len(c) < 3:
					c = '0'+c
				fig.savefig(filename+'_processed_'+c+'.png')
				fig.clf()




