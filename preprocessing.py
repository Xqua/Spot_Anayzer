#!/usr/bin/env python

import cv2
import numpy as np
import scipy.io as spio
import pims
import matplotlib.pyplot as plt


class preprocessing:

    def __init__(self, filename, matlabfile):
        self.filename = filename
        self.matlabfile = matlabfile

    def crop(self, img, mask):
        "set all pixel outside of contour to zero and resize to minimum bounding box."
        # print img.dtype, mask.dtype
        im = (img * mask) / 255
        im2 = cv2.convertScaleAbs(im)
        mask = cv2.convertScaleAbs(mask)
        # im.dtype = np.int16
        # print im.max(), im.min()
        aaa, thresh = cv2.threshold(im2, 1, 255, cv2.THRESH_BINARY)
        # thresh = cv2.convertScaleAbs(thresh)
        # print thresh.dtype, aaa, thresh.shape
        # print thresh
        # print thresh.min(), thresh.max()
        contours, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        cnt = contours[0]
        x, y, w, h = cv2.boundingRect(cnt)
        return im[y:y + h, x:x + w], mask[y:y + h, x:x + w]

    def MatLabParser(self, filepath):
        mat = spio.loadmat(filepath)
        bacterias = []
        for i in range(len(mat['Cellparsing'][0][0][0][0][0][0])):
            if len(mat['Cellparsing'][0][0][0][0][0][0][i]) > 0:
                mesh = mat['Cellparsing'][0][0][0][0][0][0][i][0][0][3]
                cnt = []
                for xyxy in mesh:
                    cnt.append([[int(round(xyxy[0], 0)), int(round(xyxy[1], 0))]])
                    cnt.append([[int(round(xyxy[2], 0)), int(round(xyxy[3], 0))]])
                bacterias.append(np.array(cnt))
        return bacterias

    def MatLab2Contour(self, mlabX, mlabY):
        cnt = []
        for i in range(len(mlabX)):
            cnt.append([[mlabX[i], mlabY[i]]])
        return cnt

    def makeMaskFromCountour(self, contour, shape):
        mask = np.zeros(shape)
        cv2.drawContours(mask, contour, -1, 255, -1)
        kernel = np.ones((7, 7), np.uint8)
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
        return mask

    def getTranslationMatrix2d(self, dx, dy):
        """
        Returns a numpy affine transformation matrix for a 2D translation of
        (dx, dy)
        """
        return np.matrix([[1, 0, dx], [0, 1, dy], [0, 0, 1]])

    def rotateImage(self, image, angle):
        """
        Rotates the given image about it's centre
        """

        image_size = (image.shape[1], image.shape[0])
        image_center = tuple(np.array(image_size) / 2)

        rot_mat = np.vstack([cv2.getRotationMatrix2D(image_center, angle, 1.0), [0, 0, 1]])
        trans_mat = np.identity(3)

        w2 = image_size[0] * 0.5
        h2 = image_size[1] * 0.5

        rot_mat_notranslate = np.matrix(rot_mat[0:2, 0:2])

        tl = (np.array([-w2, h2]) * rot_mat_notranslate).A[0]
        tr = (np.array([w2, h2]) * rot_mat_notranslate).A[0]
        bl = (np.array([-w2, -h2]) * rot_mat_notranslate).A[0]
        br = (np.array([w2, -h2]) * rot_mat_notranslate).A[0]

        x_coords = [pt[0] for pt in [tl, tr, bl, br]]
        x_pos = [x for x in x_coords if x > 0]
        x_neg = [x for x in x_coords if x < 0]

        y_coords = [pt[1] for pt in [tl, tr, bl, br]]
        y_pos = [y for y in y_coords if y > 0]
        y_neg = [y for y in y_coords if y < 0]

        right_bound = max(x_pos)
        left_bound = min(x_neg)
        top_bound = max(y_pos)
        bot_bound = min(y_neg)

        new_w = int(abs(right_bound - left_bound))
        new_h = int(abs(top_bound - bot_bound))
        new_image_size = (new_w, new_h)

        new_midx = new_w * 0.5
        new_midy = new_h * 0.5

        dx = int(new_midx - w2)
        dy = int(new_midy - h2)

        trans_mat = self.getTranslationMatrix2d(dx, dy)
        affine_mat = (np.matrix(trans_mat) * np.matrix(rot_mat))[0:2, :]
        result = cv2.warpAffine(image, affine_mat, new_image_size, flags=cv2.INTER_LINEAR)

        return result

    def GetMajorAngle(self, cnt):
        M = cv2.moments(cnt)
        u20 = M['mu20'] / M['m00']
        u02 = M['mu02'] / M['m00']
        u11 = M['mu11'] / M['m00']
        theta = 0.5 * np.arctan((2 * u11) / (u20 - u02))
        return np.degrees(theta)

    def Main(self):
        ###### PSEUDO CODE ########
        img = pims.TiffStack(self.filename)
        shape = img[0].shape
        cells = self.MatLabParser(self.matlabfile)
        c = 0
        for cnt in cells:
            result = []
            print "Extracting Bacteria %s/%s" % (c, len(cells))
            d = 0
            for im in img[:200]:
                print "Doing frame %s/%s" % (d, len(img))
                d += 1
                im = np.array(im)
                mask = np.flipud(self.makeMaskFromCountour(cnt, shape).T)
                bacteria, bacteria_mask = self.crop(im, mask)
                if d == 1:
                    fig = plt.figure()
                    plt.imshow(im)
                    plt.imshow(mask, alpha=0.5)
                    plt.show()
                    raw_input("NEXT")
                # print im.min(), im.max()
                # print bacteria.min(), bacteria.max()
                # print im.sum(), bacteria.sum()
                # print bacteria.dtype, bacteria_mask.dtype
                # print bacteria, bacteria_mask
                # im = cv2.convertScaleAbs()
                contours, hierarchy = cv2.findContours(bacteria_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
                angle = self.GetMajorAngle(contours[0])
                bacteria = self.rotateImage(bacteria, angle)
                bacteria_mask = (bacteria > 0) * 255

                bacteria, bacteria_mask = self.crop(bacteria, bacteria_mask)
                if bacteria.shape[0] > bacteria.shape[1]:
                    bacteria = bacteria.T
                result.append(bacteria)
                # plot
                # im = cv2.ellipse(im,ellipse,(0,255,0),2)
            print "Saving Bacteria"
            np.save("bacteria_%s.npy" % c, result)
            c += 1
