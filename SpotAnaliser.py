import matplotlib as mpl
import matplotlib.pyplot as plt
import pims
import numpy as np
from scipy import ndimage
from scipy.interpolate import UnivariateSpline
from scipy.stats import linregress
import scipy.io as spio
import fit
import copy
import time as libtime

mpl.rc('figure', figsize=(20, 20))
mpl.rc('image', cmap='gray', interpolation='none')


class Spot_Analysis:

    def __init__(self, filename, binning=False):
        #################################
        # Loading the data and starting the preprocessing
        #################################
        self.filename = filename
        if ".mat" in self.filename:
            self.img, self.meta = self.MatLabParser(self.filename)
            tmp = []
            for i in range(300):
                tmp.append(self.img[i].T)
            self.img = np.array(tmp)
        elif ".tif" in self.filename or ".tiff" in self.filename:
            self.img = pims.TiffStack(filename)
            self.meta = None
        self.time_pool, self.binned = self.Pool_Time(self.img, Bin=binning)
        self.max_signal = self.time_pool.max()
        print "Loaded: ", filename
        self.GMM = fit.GaussianMixture()
        #################################

        #################################
        # Definition of constants
        #################################
        self.deltax = 0.005  # Defines the delta x step for derivation etc
        self.mid = len(self.time_pool) / 2  # Define Middle pixel column
        # Defines the X pixel to use for speed analysis
        self.pix_col = range(self.mid - 1, self.mid + 2)
        # Size in nm of the pixel (here CMOS camera with 100x objective = 65nm)
        self.pix_size = 65.0
        # Time between Frames in second
        self.timeinterval = 1
        # Defined how far apart can peak be (in second)
        self.timethreshold = 6 * self.timeinterval * (1 / self.deltax)  # Classic 1 second imaging

        # If in silico data, grabs number from simulation logs.
        if self.meta:
            self.timeinterval = self.meta['sec_per_frame'][0][0]
            self.pix_size = self.meta['um_per_px'][0][0] * 1000
            self.timethreshold = 60 * self.timeinterval * (1 / self.deltax)

        # Speed Threshold, if an object moves faster than this, it is excluded
        # !
        self.speed_threshold = 40
        # Allowed variance in second for a peak to be called true when found in
        # vertical axis and checked VS time
        # self.Time_Error = 10 / self.timeinterval
        self.Time_Error = 50

        self.diameter = 850  # Diameter of the bacteria in nm

        # Minimum signal intensity to be called a peak (the signal is
        # normalized and goes from 0-1 should be set to mean(noise) +
        # 2std(noise) or learned )
        self.y_threshold = 0.06
        # Threshold value to detect the peaks > How close to zero does the
        # first derivative need to be to be considered a peak
        self.dy_threshold = 0.01
        # self.ddy_threshold = -0.1
        # Value to calculate the real length
        self.bact_lenght = len(self.time_pool[0].T)
        self.bact_lenght_nm = self.bact_lenght * self.pix_size
        # Which model to use, gaussian or flat top gaussian, if True, flat top, if False Gaussian
        self.Flat_Top = True
        # Strenght of overfit correction, the higher the more the number of parameter is penalized (1 is an equivalent Fit Vs NB of Free param)
        self.BIC_Lambda = 4
        #################################

        #################################
        # Definition of Place Holders
        #################################
        self.all_v = []
        self.All_hit = []
        self.All_time = {}
        self.All_Sigmas = []
        self.All_Intensity = []
        self.All_Spread = []
        self.Real_hit = []
        #################################

    def MatLabParser(self, filepath):
        mat = spio.loadmat(filepath)
        im = mat['im'].T
        labels = [('num_frames', 'O'),
                  ('sec_per_frame', 'O'),
                  ('exposure', 'O'),
                  ('time_step', 'O'),
                  ('diff_coeff', 'O'),
                  ('u_convection', 'O'),
                  ('um_per_px', 'O'),
                  ('sim_box_size_um', 'O'),
                  ('box_size_px', 'O'),
                  ('brightness', 'O'),
                  ('signal_background', 'O'),
                  ('num_particle', 'O'),
                  ('BD_type', 'O'),
                  ('n_dims', 'O'),
                  ('bound_condi', 'O'),
                  ('drift', 'O'),
                  ('psf_type', 'O'),
                  ('psf_sigma_um', 'O'),
                  ('renderer', 'O'),
                  ('offset', 'O'),
                  ('readout_noise', 'O'),
                  ('EMgain', 'O'),
                  ('ADCgain', 'O'),
                  ('finer_grid', 'O'),
                  ('Pe', 'O'),
                  ('num_particle_1', 'O'),
                  ('num_particle_2', 'O'),
                  ('u_convection_1', 'O'),
                  ('u_convection_2', 'O'),
                  ('density', 'O'),
                  ('bonds_per_atom', 'O'),
                  ('spring_const', 'O'),
                  ('topo', 'O'),
                  ('make_gradient', 'O'),
                  ('ic', 'O')]
        metadata = mat['o'][0][0]
        meta = {}
        for i in range(len(metadata)):
            meta[labels[i][0]] = metadata[i]
        return im, meta

    def bin_img(self, img, Sum=True, width=1, height=1):
        "Bin a 2D picture with binning of size Width and height. If Sum==True, then sums to pixels, otherwise, average them"
        res = []
        if width > 1:
            for im in img.T:
                data = im
                if Sum:
                    res.append(
                        data[:(data.size // width) * width].reshape(-1, width).sum(axis=1))
                else:
                    res.append(
                        data[:(data.size // width) * width].reshape(-1, width).mean(axis=1))
            res = np.array(res)
        else:
            res = img.T
        if height > 1:
            res2 = []
            for im in res:
                data = im
                if Sum:
                    res2.append(
                        data[:(data.size // width) * width].reshape(-1, width).sum(axis=1))
                else:
                    res2.append(
                        data[:(data.size // width) * width].reshape(-1, width).mean(axis=1))
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
        idx = (np.abs(array - value)).argmin()
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
        X = np.linspace(0, len(Y) * self.timeinterval, len(Y))
        s = linregress(X, Y)
        return s[0], s[1]

    def Signal_Process(self, Y, Gaussian=True, Gaussian_sigma=1, Log=False, MinRemove=True):
        if Gaussian:
            Y = ndimage.filters.gaussian_filter1d(
                Y, sigma=Gaussian_sigma, mode='reflect')
        if Log:
            Y = np.log(Y) / np.log(self.max_signal)
        else:
            Y = Y / float(self.max_signal)
        if MinRemove:
            Y = Y - Y.min()
        return Y

    def Peak_detect(self, Y, y_threshold=None, dy_threshold=None, point_spacing=60, smoothing_sigma=3, process=True):
        """Y is the signal to be processed.
        y_threshold defines the minimum signal amount to be a peak (filter out noise)
        dy_threshold defines the bounds in witch the derivative needs to be to call a peak
        point_spacing defines the minimum spacing between two peak too filter out multiple calls of the same peak
        smoothing_sigma defines the spread of the gaussian smoothing to preprocess the signal"""
        # Process Signal
        if process:
            Y = self.Signal_Process(Y, Gaussian_sigma=smoothing_sigma)
        if not y_threshold:
            y_threshold = self.y_threshold
        if not dy_threshold:
            dy_threshold = self.dy_threshold
        # Create X vectors
        # X = np.linspace(0, len(Y) * self.timeinterval, len(Y))
        # xs = np.arange(0, len(Y) * self.timeinterval, self.deltax)
        X = np.linspace(0, len(Y), len(Y))
        xs = np.arange(0, len(Y), self.deltax)

        # Fit Spline model to the data
        spl = UnivariateSpline(X, Y, k=5, s=0)

        # Take the Spline derivative to look for d(y) = 0 and dd(y) <= 0
        dspl = spl.derivative()  # first derivative
        ddspl = dspl.derivative()  # second derivative

        # generate Smooth data from the curve spline fitting
        y = spl(xs)
        dy = dspl(xs)
        ddy = ddspl(xs)

        # Find the roots, dy == 0 and ddy <= 0
        roots = []
        for i in range(len(dy)):
            if (-dy_threshold < dy[i] < dy_threshold) and y[i] > y_threshold and ddy[i] <= 0:
                # print y[i]
                roots.append(i)

        # Filter out the points that are next to each other
        roots2 = []
        for i in range(len(roots)):
            if i < len(roots) - 1:
                if not (roots[i + 1] - roots[i]) < point_spacing:
                    roots2.append(roots[i])
            else:
                roots2.append(roots[i])
        roots = roots2
        return roots, X, xs, Y, y, dspl(xs), ddspl(xs)

    def Filter_Neighbor_Hit(self, roots_all, xs):
        hitlist = []
        roots_ok = []
        for r in roots_all[1]:
            for r0 in roots_all[0]:
                for r2 in roots_all[2]:
                    if abs(r - r2) < self.timethreshold and abs(r - r0) < self.timethreshold and abs(r2 - r0) < 2 * self.timethreshold:
                        if r0 < r < r2 or r2 < r < r0:
                            roots_ok.append((r0, 'b'))
                            roots_ok.append((r, 'g'))
                            roots_ok.append((r2, 'r'))
                        else:
                            roots_ok.append((r0, 'k'))
                            roots_ok.append((r, 'k'))
                            roots_ok.append((r2, 'k'))
                        t = np.mean(
                            (abs(xs[r] - xs[r0]), abs(xs[r] - xs[r2]), abs(xs[r0] - xs[r2]) / 2))
                        v = self.pix_size / t
                        if v < self.speed_threshold:
                            self.all_v.append(v)
                        hitlist.append(int(round(xs[r], 0)))
        return roots_ok, hitlist

    def Fit(self, X, xs, Y, roots, backfit):
        rr = np.array([xs[i] for i in roots])
        F = self.GMM.BackFit(X, Y, rr, fusion=self.Flat_Top, backfit=backfit)
        BIC = self.GMM.BIC(X, Y, [[i.x, 0] for i in F], fusion=self.Flat_Top, Lambda=self.BIC_Lambda)
        b_i, b = np.argmin(BIC), BIC[np.argmin(BIC)]
        print 'BICs : ', BIC
        print "Best BIC: ", b
        model = F[b_i]
        if self.Flat_Top:
            roots, sigma, amplitude, spread = np.reshape(copy.copy(model['x']), (len(model['x']) / 4, 4)).T
            # print roots, spread
            roots += (spread / 2)
            # print roots, model['x']
        else:
            roots, sigma, amplitude = np.reshape(model['x'], (len(model['x']) / 3, 3)).T
            spread = [0] * len(roots)
        roots = np.searchsorted(xs, roots) - 1
        return roots, sigma, amplitude, spread, model

    def Main(self, makeplot=False):
        t0 = libtime.time()
        time_hit = []
        color = {self.pix_col[0]: 'b', self.pix_col[1]: 'g', self.pix_col[2]: 'r'}
        print "Bacteria Lenght:", len(self.time_pool[0].T)
        ddys = {}

        for line in range(len(self.time_pool[0].T)):
            print "line is :", line
            roots_all = []
            # sigma_all = []
            # intensity_all = []
            # spread_all = []

            if makeplot:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_ylim(ymax=1)
            roots_max = []
            Ys = []
            for tp in self.pix_col:
                Y = self.time_pool[tp].T[line]
                roots, X, xs, Y, y, dy, ddy = self.Peak_detect(Y)
                # print Y
                # print roots
                Ys.append(Y)
                if len(roots) > len(roots_max):
                    roots_max = roots
            if len(roots_max) > 0:
                Y = np.mean(Ys, axis=0)
                # print len(Y), Y
                roots_start, sigma, amplitude, spread, model = self.Fit(X, xs, Y, roots_max, True)
            else:
                roots_start = []
            print "Root Bootsart:", len(roots_start)
            for tp in self.pix_col:
                Y = self.time_pool[tp].T[line]
                roots, X, xs, Y, y, dy, ddy = self.Peak_detect(Y)
                # print 'roots_1', len(roots)
                if len(roots_start) > 0:
                    roots, sigma, amplitude, spread, model = self.Fit(X, xs, Y, roots_start, False)
                    print 'roots_2', len(roots)
                    [self.All_Sigmas.append(i) for i in sigma]
                    [self.All_Intensity.append(i) for i in amplitude]
                    [self.All_Spread.append(i) for i in spread]
                    if self.pix_col.index(tp) == 1:
                        ddys[line] = ddy

                    if makeplot:
                        self.GMM.Plot(X, Y, model['x'], fusion=self.Flat_Top, ax=ax, color=color[tp])
                else:
                    roots = []
                if makeplot:
                    ax.plot(X, Y, label='pix%s' % tp, color=color[tp])
                    # ax.plot(xs, y, label='pix%s' % tp, color=color[tp])
                    # ax.plot(xs, dy, 'orange')
                    # ax.plot(xs, ddy, 'c')
                    # ax.plot(xs, [0] * len(xs), 'y')
                roots_all.append(roots)

            print 'roots_all', len(roots_all)
            roots_ok, hitlist = self.Filter_Neighbor_Hit(roots_all, xs)
            print "hitlist: ", hitlist
            print 'roots_ok', len(roots_ok)
            self.All_hit.append(np.unique(hitlist))
            # Plot vertical line on the dectected peaks
            for i in roots_ok:
                if makeplot:
                    plt.plot((xs[i[0]], xs[i[0]]), (0, 1), i[1])
                if i[1] == 'g' or i[1] == 'k':
                    time_hit.append(self.find_nearest(X, xs[i[0]]))

            if makeplot:
                plt.legend()
                fig.savefig(self.filename.split('.')[0] + '_time_peak_detect' + str(line) + '.png')
                fig.clf()
                plt.clf()
                del fig

            # raw_input("NEXT")
        time_hit = np.unique(time_hit)
        time_hit = [np.where(X == i)[0][0] for i in time_hit]
        self.speed = np.mean(self.all_v)
        if np.isnan(self.speed):
            self.speed = 0
        print "mean Speed", self.speed
        for time in time_hit:
            # print "time is :",time
            Y = self.Signal_Process(self.time_pool[self.mid][time], Gaussian=False)
            # print Y
            roots, X, xs, Y, y, dy, ddy = self.Peak_detect(Y, process=False)
            # print roots
            # print Y
            # self.Flat_Top = False
            # roots, sigma, amplitude, spread, model = self.Fit(X, xs, Y, roots, True)
            # self.Flat_Top = True
            self.All_time[time] = [int(round((self.find_nearest(X, xs[i])) - 1, 0)) for i in roots]

            # Plot vertical line on the dectected peaks
            if makeplot:
                fig = plt.figure()
                plt.plot(X, Y)
                plt.plot(xs, y)
                plt.plot(xs, dy, 'orange')
                plt.plot(xs, ddy, 'c')
                plt.plot(xs, [0] * len(xs), 'y')
                for i in roots:
                    plt.plot((xs[i], xs[i]), (0, 1))
                plt.legend()
                fig.savefig(self.filename.split('.')[0] + '_vertical_peak_detect' + str(time) + '.png')
                fig.clf()
                plt.clf()
                del fig

        self.Real_hit = []
        self.gauss_sigma = []
        for t in self.All_time.keys():
            for line in self.All_time[t]:
                test = False
                for l in range(line - 1, line + 2):
                    try:
                        n = self.find_nearest(self.All_hit[line], t)
                        if abs(n - t) < self.Time_Error and n >= 0:
                            test = True
                    except:
                        pass
                if test:
                    self.Real_hit.append((line + 1, t))
                    # self.gauss_sigma.append(ddys[line + 1][n])
        if self.Real_hit:
            Tn = {}
            lines, times = np.array(self.Real_hit)[:, 0], np.array(self.Real_hit)[:, 1]
            for i in range(len(times)):
                if times[i] in Tn:
                    Tn[times[i]].append(lines[i])
                else:
                    Tn[times[i]] = [lines[i]]
            T = np.unique(times)
            res = np.zeros(len(T))

            for t in range(len(T)):
                c = 0
                for l in Tn[T[t]]:
                    if t < len(T) - 2:
                        # if l not in Tn[T[t + 2]] and l not in Tn[T[t + 1]]:
                        if (l not in Tn[T[t + 2]] and l not in Tn[T[t + 1]]) and (l + 1 not in Tn[T[t + 2]] and l + 1 not in Tn[T[t + 1]]) and (l - 1 not in Tn[T[t + 2]] and l - 1 not in Tn[T[t + 1]]):
                            c += 1
                    else:
                        c += 1
                res[t] = c

            print T, res

            # self.revolution_window = int((self.diameter * np.pi) / self.speed)
            # windowed = []
            # for i in range(len(res) - self.revolution_window):
            #     windowed.append(np.sum(res[i:i + self.revolution_window]))
            self.meanMREB = np.sum(res)
            self.stdMREB = 0
        else:
            self.meanMREB, self.stdMREB = 0, 0
        print "Number of Events per revolution window:"
        print "Mean:", self.meanMREB, "Std:", self.stdMREB
        print "Gaussian width:"
        # print "Mean:", (np.sqrt(np.abs(np.mean(1.0 / np.array(self.gauss_sigma))))), "std:", (np.sqrt(np.abs(np.std(1.0 / np.array(self.gauss_sigma)))))
        print "Putative filament Lenght:"
        print (np.sqrt(np.abs(np.mean(1.0 / np.array(self.gauss_sigma)))) * self.speed) - 240
        print "Time Elapsed for the bacteria in second: ", libtime.time() - t0

        if makeplot:
            # fig = plt.figure()
            # plt.plot(range(len(windowed)), windowed)
            # plt.title('Number of MReB throught the sliding window, Mean=%s, Std=%s' % (
            #     self.meanMREB, self.stdMREB))
            # fig.savefig(self.filename.split('.')[0] + 'MReB_Detection.png')
            # fig.clf()
            for i in range(len(self.binned)):
                fig = plt.figure()
                plt.imshow(self.binned[i])
                tt = np.where(times == i)[0]
                plt.plot(0, 0, 'k.')
                for t in tt:
                    plt.plot(self.mid, lines[t], 'ro')
                c = str(i)
                while len(c) < 3:
                    c = '0' + c
                fig.savefig(self.filename + '_processed_' + c + '.png')
                fig.clf()
                plt.clf()
                del fig
