##########################################
#                                        #
#         ROBOTIC SCIENCE AND SYSTEMS    #
#         Spring 2015                    #
#         Prof: Rob Platt                #
#                                        #
#         Author: Kirby Powell           #
#         Due: 04-29-2015                #
#                                        #
##########################################

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from math import sqrt

def plot_cov_ellipse(cov, pos, volume=.5, ax=None, fc='none', ec=[0,0,0], a=1, lw=2):
    """
    Plots an ellipse enclosing *volume* based on the specified covariance
    matrix (*cov*) and location (*pos*). Additional keyword arguments are passed on to the
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        volume : The volume inside the ellipse; defaults to 0.5
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.

    Code found at:
    http://www.nhsilbert.net/source/2014/06/bivariate-normal-ellipse-plotting-in-python/
    """

    import numpy as np
    from scipy.stats import chi2
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    kwrg = {'facecolor':fc, 'edgecolor':ec, 'alpha':a, 'linewidth':lw}

    # Width and height are "full" widths, not radius
    width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwrg)

    ax.add_artist(ellip)


def extractXY(answers):
	"""Takes a list of kalman-filtered values and extracts the first two values from each
	value, expecting those two to be x and y values"""

	xy_answers = []
	for i in range(0, len(answers)):
		xy_answers.append((answers[i].item(0), answers[i].item(1)))

	return xy_answers


def separateXY(xy_answers):
	"""returns a list of x and y values from a given list of x,y tuples"""
	x = []
	y = []

	for i in range(0, len(xy_answers)):
		x.append(xy_answers[i][0])
		y.append(xy_answers[i][1])

	return x,y

def extractCovariance(listofSigmas):
	"""returns a list of 2x2 covariance matrices from a given list 
	NxN covariance matrices"""

	listOfCovs = []

	for sigma in listofSigmas:
		listOfCovs.append(np.matrix([[sigma.item((0,0)), sigma.item((0,1))], [sigma.item((1,0)), sigma.item((1,1))]]))

	return listOfCovs

def partialExtendedKalman(A, Q, init_mean, init_cov, observes):
	"""The Kalman filter method. Takes all expected matrices,
	along with a list of observations (expected to be a list
	of matrices), and returns the Kalman filtered values in a list."""

	mu = init_mean
	sigma = init_cov

	# since the initial mean is the first value, store it in the list of mus
	listOfMuPrime = [mu]
	listofSigPrime = [sigma]

	# loop through the observations and record filtered values
	for i in range(1,len(observes)):
		# PROCESS UPDATE
		muPrime = A * mu #"""THIS MUST CHANGE FOR EKF"""
		sigPrime = (A * sigma * A.transpose()) + Q

		mu = muPrime
		sigma = sigPrime

		# save the filtered value & repeat
		listOfMuPrime.append(muPrime)
		listofSigPrime.append(sigPrime)

	return listOfMuPrime, listofSigPrime

def fullExtendedKalman(A, R, Q, init_mean, init_cov, observes):
	"""The Kalman filter method. Takes all expected matrices,
	along with a list of observations (expected to be a list
	of matrices), and returns the Kalman filtered values in a list."""

	# set up mu & sigma
	mu = init_mean
	sigma = init_cov

	# since the initial mean is the first value, store it in the list of mus
	listOfMu = [mu]
	listofSigPrime = [sigma]

	# loop through the observations and record filtered values
	for i in range(1,len(observes)):
		# PROCESS UPDATE
		muPrime = A * mu #"""THIS MUST CHANGE FOR EKF"""
		sigPrime = (A * sigma * A.transpose()) + Q

		C = np.matrix([[(((muPrime.item(1)**2 + muPrime.item(0)**2)**(-1/2))*muPrime.item(0)), muPrime.item(1)*(muPrime.item(1)**2 + muPrime.item(0)**2)**(-1/2), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
		H = np.matrix([[sqrt(muPrime.item(1)**2 + muPrime.item(0)**2)],[muPrime.item(2)], [muPrime.item(3)]])

		# MEASUREMENT UPDATE
		mu = muPrime + sigPrime*C.transpose() * inv(R + C * sigPrime * C.transpose()) * (observes[i] - H)
		sigma = sigPrime - sigPrime*C.transpose() * inv(R + C * sigPrime * C.transpose()) * C * sigPrime

		# save the filtered value & repeat
		listOfMu.append(mu)
		listofSigPrime.append(sigPrime)

	return listOfMu, listofSigPrime

def main():
	"""Main method. Calls all above methods directly or indirectly, 
	and plots all applicable graphs."""

	observations = [np.matrix('3.1623; 0; .5'), np.matrix('2.6926; 0; .5'), np.matrix('2.2361; 0; .5'), np.matrix('1.8028; 0; .5'),
				np.matrix('1.4142; 0; .5'), np.matrix('1.1180; 0; .5'), np.matrix('1.0000; 0; .5'), np.matrix('1.1180; 0; .5'),
				np.matrix('1.4142; 0; .5'), np.matrix('1.8028; 0; .5'), np.matrix('2.2361; 0; .5'), np.matrix('2.6926; 0; .5'), np.matrix('3.1623; 0; .5')]
	priorMean = np.matrix('2; -3; 0; .5')

	A = np.matrix('1 0 1 0; 0 1 0 1; 0 0 0.9 0; 0 0 0 0.9')
	prirorCovariance = np.matrix('1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1')
	Q = np.matrix('.1 0 0 0; 0 .1 0 0; 0 0 .1 0; 0 0 0 .1')
	R = np.matrix('.1 0 0 ; 0 .1 0 ; 0 0 .1 ')

	partialAnswersMu, partialAnswersSigma = partialExtendedKalman(A, Q, priorMean, prirorCovariance, observations)

	fullAnswersMu, fullAnswersSigma = fullExtendedKalman(A, R, Q, priorMean, prirorCovariance, observations)

	partialMuXY = extractXY(partialAnswersMu)
	partialMuX, partialMuY = separateXY(partialMuXY)

	fullMuXY = extractXY(fullAnswersMu)
	fullMuX, fullMuY = separateXY(fullMuXY)

	partialCovariances = extractCovariance(partialAnswersSigma)
	fullCovariances = extractCovariance(fullAnswersSigma)

	plt.plot(partialMuX, partialMuY, "bx-")
	plt.title('V vs W\nProcess Update Only')
	plt.ylabel('W Axis')
	plt.xlabel('V Axis')

	for i in range(len(partialCovariances)):
		plot_cov_ellipse(partialCovariances[i] ,partialMuXY[i], .00001)

	plt.show()

	plt.plot(fullMuX, fullMuY, "bx-")
	plt.title('V vs W\nProcess and Measurement Updates')
	plt.ylabel('W Axis')
	plt.xlabel('V Axis')

	for i in range(len(fullCovariances)):
		plot_cov_ellipse(fullCovariances[i] ,fullMuXY[i], .01)

	plt.show()

if __name__ == '__main__':
	main()