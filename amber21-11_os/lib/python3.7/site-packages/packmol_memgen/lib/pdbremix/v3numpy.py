# encoding: utf-8

__doc__ = """

3D vector geometry library based on numpy.

The vectors and transform matrices are subclassd from numpy.array. 
All operations are accessed through functions, allowing easy
switching with other libraries, such as v3array.
"""

import math
import random

import numpy as np


def vector(*args):
  """
  Creates a new vector

  Args:
    None: returns zero vector
    v (vector): returns copy of v
    x, y, z: returns (x, y, z)
  """
  n_arg = len(args)
  if n_arg == 0:
    return np.zeros(3, dtype=np.float)
  if n_arg == 1:
    data = args[0]
    if len(data) == 3:
      return np.array(data, copy=True)
    raise TypeError('vector() with 1 argument must have 3 elements')
  if n_arg == 3:
    return np.array(args, dtype=np.float, copy=True)
  else:
    raise TypeError('vector() takes 0, 1 or 3 arguments')


def set_vector(*args):
  """
  Changes values of vector in place

  Args:
    v (vector), w (vector): copies values of w into v
    v (vector), x, y, z: copiex x, y, z into v
  """
  vector = args[0]
  if len(args) == 2:
    vector[:] = args[1]
  elif len(args) == 4:
    vector[:] = args[1:4]


def mag(vector):
  """
  Returns the magnitude of a vector
  """
  return np.sqrt(np.dot(vector, vector))


def mag2(vector):
  """
  Returns square of the magnitude of a vector
  """
  return np.dot(vector, vector)


def scale(vector, s):
  """
  Returns vector that has been scaled by s.
  """
  return  s*vector


dot = np.dot


def norm(vector):
  """
  Returns vector normalized to magnitude=1.0
  """
  return vector/mag(vector)


cross = np.cross


radians = np.radians


degrees = np.degrees


def identity():
  """
  Returns the identity transform.
  """
  return np.eye(4)


def matrix_elem(matrix, i, j, val=None):
  """
  Reads/writes the elements of an affine transform.

  1. 3x3 rotational component;
      matrix_elem(m, i, j) for i=0..3, j=0..3
  2. 3x1 translational component:
      matrix_elem(m, 3, i) for i=0..3)
  """
  if val is None:
    return matrix[j,i]
  else:
    matrix[j,i] = val


def transform(matrix, vector):
  """
  Returns vector of applying the transform in matrix to v.
  """
  # return np.dot(matrix[:3,:3], vector) + matrix[:3,3]  
  return np.dot(matrix, np.append(vector,[1],0))[:3]


def left_inverse(matrix):
  """
  Returns the left inverse of m.

  Example:
    combine(left_inverse(m), m) == identity().
  """
  inverse = identity()
  r = matrix[:3,:3].transpose()
  inverse[:3,:3] = r
  inverse[:3,3] = -np.dot(r, matrix[:3,3])
  return inverse


def rotation(axis, theta):
  """
  Returns transform that rotate a vector at the origin around axis.
  """
  # from http://stackoverflow.com/a/6802723
  m = identity()
  a = np.cos(theta/2)
  b, c, d = norm(axis) * np.sin(theta/2)
  m[0,:3] = [a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)    ]
  m[1,:3] = [2*(b*c+a*d),     a*a+c*c-b*b-d*d, 2*(c*d-a*b)    ]
  m[2,:3] = [2*(b*d-a*c),     2*(c*d+a*b),     a*a+d*d-b*b-c*c]
  return m


def translation(displacement):
  """
  Returns transform that translates a vector.
  """
  m = identity()
  m[:3,3] = displacement
  return m


def scaling_matrix(s0, s1, s2):
  """
  Returns a scaling matrix.
  """
  m = identity()
  matrix_elem(m, 0, 0, s0)
  matrix_elem(m, 1, 1, s1)
  matrix_elem(m, 2, 2, s2)
  return m


def combine(m1, m2):
  """
  Returns transform that combines two other transforms.
  """
  return np.dot(m1, m2)


def is_similar_mag(a, b, small=1E-5):
  """
  Evaluates similar magnitudes to within small.
  """
  return abs(abs(a)-abs(b)) <= small


def is_similar_matrix(m1, m2, small=1E-5):
  """
  Evaluates similar matrixes through matrix components.
  """
  iter1 = np.ndenumerate(m1)
  iter2 = np.ndenumerate(m2)
  for (i1, val1), (i2, val2) in zip(iter1, iter2):
    if not is_similar_mag(val1, val2, small):
      return False
  return True


is_similar_vector = is_similar_matrix


# common to v3list and v3numpy


def parallel(v, axis):
  """
  Returns component of v parallel to axis.
  """
  l = mag(axis)
  if is_similar_mag(l, 0):
    return v
  else:
    return scale(axis, dot(v, axis)/l/l) 


def perpendicular(v, axis):
  """
  Returns component of v that is perpendicular to axis.
  """
  return v - parallel(v, axis)


def normalize_angle(angle):
  """
  Returns angle in radians that is [-pi, pi]
  """
  while abs(angle) > math.pi:
    if angle > math.pi:
      angle -= math.pi*2
    if angle < -math.pi:
      angle += 2*math.pi
  if is_similar_mag(abs(angle + math.pi), 0):
    angle = math.pi
  return angle


def vec_angle(a, b):
  """
  Returns angle in radians between a and b. 
  """ 
  a_len = mag(a)
  b_len = mag(b)
  if is_similar_mag(a_len, 0) or is_similar_mag(b_len, 0):
    return 0.0
  c = dot(a, b) / a_len / b_len
  if c >=  1.0:
    return 0.0
  elif c <= -1.0:
    return math.pi
  else:
    return math.acos(c)  


def vec_dihedral(a, axis, c):
  """
  Returns dihedral angle between a and c, along the axis.
  """ 
  ap = perpendicular(a, axis)
  cp = perpendicular(c, axis)
  angle = vec_angle(ap, cp)
  if dot(cross(ap, cp), axis) > 0:
    angle = -angle
  return angle


def dihedral(p1, p2, p3, p4):
  """
  Returns dihedral angle defined by the four positions.
  """ 
  return vec_dihedral(p1-p2, p2-p3, p4-p3)


def distance(p1, p2):
  """
  Returns distance between the two points
  """ 
  return mag(p1 - p2)


def rotation_at_center(axis, theta, center):
  """
  Returns a rotation around the center.
  """ 
  t = translation(-center)
  r = rotation(axis, theta)
  t_inv = translation(center)
  return combine(t_inv, combine(r, t))


def get_center(crds):
  """
  Returns the geometric center of a bunch of positions.
  """
  center = vector()
  for crd in crds:
    center += crd
  return scale(center, 1.0/float(len(crds)))


def get_width(crds):
  """
  Returns the maximum width between any two crds in the group.
  """
  center = get_center(crds)
  max_diff = 0
  for crd in crds:
    diff = v3.distance(crd, center)
    if diff > max_diff:
      max_diff = diff
  return 2*max_diff


def random_mag(): 
  """
  Returns a random positive number from [0, 90] for testing.
  """
  return random.uniform(0, 90)


def random_real(): 
  """
  Returns a random real +/- from [-90, 90] for testing.
  """
  return random.uniform(-90, 90)


def random_vector():
  """
  Returns a random vector for testing.
  """
  return vector(random_real(), random_real(), random_real())


def random_rotation():
  """
  Returns a random rotational matrix for testing.
  """
  return rotation(random_vector(), radians(random_real()))


def random_matrix():
  """
  Returns a random transformation matrix for testing.
  """
  return combine(random_rotation(), translation(random_vector()))




