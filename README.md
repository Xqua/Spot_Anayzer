# Spot_Anayzer

USAGE: run.py PATH_TO_IMAGE_FILES TXT_IMAGE_FILE_LIST
	to run the example:
	python run.py example/images list.txt

The file SpotAnalyser.py contains in the __init__ function lots of parameter that need to be changed if need be:

		self.deltax = 0.1 						#Defines the delta x step for derivation etc
		self.mid = len(self.time_pool)/2 			#Define Middle pixel column
		self.pix_col = range(self.mid-1,self.mid+2)		#Defines the X pixel to use for speed analysis
		self.timeinterval = 2 					# Time between Frames in second 
		self.speed_threshold = 40 				#Speed Threshold, if an object moves faster than this, it is excluded !
		self.timethreshold = 4 *self.timeinterval*(1/self.deltax) #Defined how far appart can peak be (in second)
		self.Time_Error = 6 / self.timeinterval #Allowed variance in second for a peak to be called true when found in vertical axis and checked VS time  
		self.diameter = 850 					#Diameter of the bacteria in nm
		self.pix_size = 65 						#Size in nm of the pixel (here CMOS camera with 100x objective = 65nm)
