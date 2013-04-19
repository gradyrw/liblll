# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Based on the implementation of liblll from https://github.com/kutio/liblll
# Modified by trdean to use numpy.  Got rid of use of Fraction...computing
# gcd on each call was taking most of time to compute LLL (now using floats)

import numpy
import math

# gram-schmidt algorithm
def gram_schmidt(g, m, mu, B):

  row = len(m)

  for i in xrange(row):
    # bi* = bi
    b_i = g[:,i]
    b_i_star = b_i
    m[:, i] = b_i_star

    for j in xrange(i):
      # u[i][j] = (bi, bj*)/Bj
      b_j_star = m[:, j]
      b_i = g[:, i]
      B[j] = numpy.vdot(b_j_star, b_j_star)
      mu[i,j] = numpy.vdot(b_i, b_j_star)[0,0]/ float(B[j,0])
      b_i_star = b_i_star - b_j_star*mu[i,j]
      m[:, i] = b_i_star

    b_i_star = m[:, i]
    # B[i] = (bi*, bi*)
    B[i] = numpy.vdot(b_i_star, b_i_star)
 

# Performs main computation of step 1 of the LLL algorithm (achieves 1.18)
def reduce(g, mu, k, l):
  if math.fabs(mu[k,l]) > 0.5:
    r = int(round(mu[k,l]))
    b_k = g[:, k]
    b_l = g[:, l]
    # bk = bk - r*bl
    g[:, k] = b_k - b_l*r

    for j in xrange(l):
      mu[k,j] = mu[k,j] - r*mu[l,j]

    mu[k,l] = mu[k,l] - r


# Main LLL Reduction
def lll_reduction(n, lc=0.75):

  dim = n.shape[0]

  m = numpy.matrix(numpy.zeros(shape=(dim,dim)),int).astype(object)
  mu = numpy.matrix(numpy.zeros(shape=(dim,dim)))
  b = n
  B = numpy.matrix(numpy.zeros(shape=(dim,1)),int).astype(object)

  gram_schmidt(b, m, mu, B)

  #k starts at 1 instead of 2 because of zero indexing
  k = 1

  while 1:
    # Performs step 1 of the LLL algorithm
    reduce(b, mu, k, k-1)

    # Performs step 2 case 1, where case 1 
    #=> condition 1.20 from LLL paper
    if (B[k,0] < lc*B[k-1,0] - mu[k,k-1]**2 * B[k-1,0] and k > 0):

      #First exchange b_k and b_(k-1)
      tmp = numpy.matrix(b[:, k])
      b[:,k] = b[:, k-1]
      b[:,k-1] = tmp
      #all other values of b remain the same

      # Store mu[k,k-1], it will be needed later
      u = mu[k,k-1]
      
      #Next adjust the set B (referred to as C* in LLL)
      B_temp = B[k,0]
      C = B[k,0] + mu[k,k-1]**2 * B[k-1,0]
      mu[k,k-1] = mu[k,k-1]*B[k-1,0]/(1.0*C)
      B[k,0] = B[k-1,0] * B[k,0]/C
      B[k-1,0] = C
      #All other B[i] remains the same

      #Adjust the set mu next
      #for  i > k
      for i in xrange(k+1, dim):
        temp = mu[i,k-1]
        print C[k-1,0]
        mu[i,k-1] = mu[i,k-1]*mu[k,k-1] + mu[i,k]*B_temp/(1.0*B[k-1,0])
        #recall that u is the original mu[i,k]
        mu[i,k] = temp - mu[i,k]*u

      # for j = 0 .. k-1
      for j in xrange(k-1):
        temp = mu[k-1,j]
        mu[k-1,j] = mu[k,j]
        mu[k,j] = temp
      
      #All other values of mu remain the same
      k = k - 1
      print k

    else:
      l = k-1
      while l >= 0:
        if math.fabs(mu[k,l]) > 0.5: 
          reduce(b, mu, k, l)
          l = k-1
        else:
          l = l - 1
      k = k + 1
      
    print k
    if (k == dim):
      return b

# Checks if in form of LLL reduction
def islll(n, lc=0.75):

  dim = n.shape[0]

  m = numpy.matrix(numpy.zeros(shape=(dim,dim)),int).astype(object)
  mu = numpy.matrix(numpy.zeros(shape=(dim,dim)))
  B = numpy.matrix(numpy.zeros(shape=(dim,1)),int).astype(object)

  gram_schmidt(n, m, mu, B)

  for i in xrange(dim):
    for j in xrange(i):
      if math.fabs(mu[i,j]) > 0.5:
        print "Hello"
        return False

  for k in xrange(1, dim):  
    if B[k] < (lc*B[k-1,0] - B[k-1,0]*mu[k,k-1]*mu[k,k-1]):
      return False

  return True
