# encoding: utf-8

__doc__ = """

3D vector geometry library in pure Python.

The vectors and transform matrices are subclassd from array.array. 
Most operations are accessed through functions, allowing easy
switching with other libraries, such as v3numpy.
"""

import math
import random
from array import array


class Vec3(array):
  """
  Represents a 3-dimensional vector/coordinate.

  Arithmetic operators are overloaded for artihmetic vector
  operations. Scaling and other non-arithmetic operators are
  carried out with functions.
  """
  
  def __new__(cls, x=0, y=0, z=0):
    return array.__new__(cls, 'd', (x,y,z))

  def __repr__(self):
    return "Vec3(%f, %f, %f)" % (self[0], self[1], self[2])

  def __add__(self, rhs):
    return Vec3(self[0]+rhs[0], self[1]+rhs[1], self[2]+rhs[2])

  def __iadd__(self, rhs):
    self[0] += rhs[0]
    self[1] += rhs[1]
    self[2] += rhs[2]
    return self

  def __isub__(self, rhs):
    self[0] -= rhs[0]
    self[1] -= rhs[1]
    self[2] -= rhs[2]
    return self

  def __sub__(self, rhs):
    return Vec3(self[0]-rhs[0], self[1]-rhs[1], self[2]-rhs[2])

  def __neg__(self):
    return Vec3(-self[0], -self[1], -self[2])

  def __pos__(self):
    return self

  def __mul__(self, s):
    return Vec3(self[0]*s, self[1]*s, self[2]*s)

  def __rmul__(self, s):
    return Vec3(self[0]*s, self[1]*s, self[2]*s)

  def __deepcopy__(self, memo):
    cls = self.__class__
    return array.__new__(cls, 'd', (self[0],self[1],self[2]))
    

def vector(*args):
  """
  Creates a new vector

  Args:
    None: returns zero vector
    v (vector): returns copy of v
    x, y, z: returns (x, y, z)
  """
  if len(args) == 0:
    return Vec3(0, 0, 0)
  elif len(args) == 1:
    v = args[0]
    return Vec3(v[0], v[1], v[2])
  elif len(args) == 3:
    return Vec3(*args)
  else:
    raise TypeError('vector() can take 0, 1, 3 arguments only')


def set_vector(*args):
  """
  Changes components of vector in place

  Args:
    v (vector), w (vector): copies values of w into v
    v (vector), x, y, z: copiex x, y, z into v
  """
  v = args[0]
  if len(args) == 2:
    w = args[1]
    v[:] = w
  elif len(args) == 4:
    v[0], v[1], v[2] = args[1:4]


def mag(v):
  """
  Returns the magnitude of v.
  """
  x, y, z = v
  return math.sqrt(x*x + y*y + z*z)


def mag2(v):
  """
  Returns square of the magnitude of v.
  """
  x, y, z = v
  return x*x + y*y + z*z


def scale(v, s):
  """
  Returns vector that is v scaled (multiplied) by s
  """
  return vector(s*v[0], s*v[1], s*v[2])


def dot(v1, v2):
  """
  Returns the dot product of two vectors
  """
  x1, y1, z1 = v1
  x2, y2, z2 = v2
  return x1*x2 + y1*y2 + z1*z2


def norm(v):
  """
  Returns scaled v, normalized to magnitude 1.0.
  """
  return scale(v, 1.0/mag(v))


def cross(v1, v2):
  """
  Returns the cross product of two vectors. 
  """
  x1, y1, z1 = v1
  x2, y2, z2 = v2
  return vector(
      y1*z2 - z1*y2,
      z1*x2 - x1*z2,
      x1*y2 - y1*x2)


def radians(degrees):
  """
  Converts degrees to radians, which is used in math functions.
  """
  return math.pi / 180.0 * degrees


def degrees(radians):
  """
  Converts radians to degrees, better for reporting.
  """
  return 180.0 / math.pi*degrees * radians


class Matrix3d(array):
  """
  Represents affine transforms in 3D space. 

  To be accessed by the auxilary function matrix_elem(). Matrix3d
  subclasses an array of 12 floats. This is subclassed mainly to
  provide a more informative __repr_function.  
  """

  def __new__(cls, matrix_array=(1,0,0,0,1,0,0,0,1,0,0,0)):
    return array.__new__(cls, 'd', matrix_array)

  def __repr__(self):
    def str3(x, y, z): 
      return "% 10.5f, % 10.5f, % 10.5f" % (x, y, z)
    s = "Matrix3d((" + str3(*self[:3]) + ",\n"
    s += "          " + str3(*self[3:6]) + ",\n"
    s += "          " + str3(*self[6:9]) + ",\n"
    s += "          #-----------------------------------\n"
    s += "          " + str3(*self[9:12]) + "))"
    return s
    


def identity():
  """
  Returns the identity transform.
  """
  return Matrix3d()


def matrix_elem(matrix, i, j, val=None):
  """
  Reads/writes the elements of an affine transform.

  1. 3x3 rotational component;
      matrix_elem(m, i, j) for i=0..3, j=0..3
  2. 3x1 translational component:
      matrix_elem(m, 3, i) for i=0..3)
  """
  k = i*3 + j
  if val is None:
    return matrix[k]
  else:
    matrix[k] = val
    return val


def transform(matrix, v):
  """
  Returns vector of applying the transform in matrix to v.
  """
  v_x, v_y, v_z = v
  x = matrix_elem(matrix, 0, 0) * v_x + \
      matrix_elem(matrix, 1, 0) * v_y + \
      matrix_elem(matrix, 2, 0) * v_z + \
      matrix_elem(matrix, 3, 0)
  y = matrix_elem(matrix, 0, 1) * v_x + \
      matrix_elem(matrix, 1, 1) * v_y + \
      matrix_elem(matrix, 2, 1) * v_z + \
      matrix_elem(matrix, 3, 1)
  z = matrix_elem(matrix, 0, 2) * v_x + \
      matrix_elem(matrix, 1, 2) * v_y + \
      matrix_elem(matrix, 2, 2) * v_z + \
      matrix_elem(matrix, 3, 2)
  return vector(x, y, z)


def left_inverse(m):
  """
  Returns the left inverse of m.

  Example:
    combine(left_inverse(m), m) == identity().
  """
  rot_t = identity()
  for i in range(0, 3):
    for j in range(0, 3):
      val = matrix_elem(m,j,i)
      matrix_elem(rot_t,i,j,val)
  t = vector(matrix_elem(m,3,0), matrix_elem(m,3,1), matrix_elem(m,3,2))
  t_inv = transform(rot_t, t)
  for i in range(0, 3):
    matrix_elem(rot_t, 3, i, -t_inv[i])
  return rot_t


def rotation(axis, theta):
  """
  Returns transform that rotate a vector at the origin around axis.
  """
  m = identity()
  c = math.cos(theta)
  s = math.sin(theta)
  t = 1.0 - c
  x, y, z = norm(axis)
  matrix_elem(m, 0, 0, t*x*x       + c)
  matrix_elem(m, 0, 1, t*x*y + z*s    )
  matrix_elem(m, 0, 2, t*x*z - y*s    )
  matrix_elem(m, 1, 0, t*y*x - z*s    )
  matrix_elem(m, 1, 1, t*y*y       + c)
  matrix_elem(m, 1, 2, t*y*z + x*s    )
  matrix_elem(m, 2, 0, t*z*x + y*s    )
  matrix_elem(m, 2, 1, t*z*y - x*s    )
  matrix_elem(m, 2, 2, t*z*z       + c)
  return m


def translation(displacement):
  """
  Returns transform that translates by displacement.
  """
  m = identity()
  for i in range(3):
    matrix_elem(m, 3, i, displacement[i])
  return m


def combine(a, b):
  """
  Returns transform that combines two transforms.
  """
  c = identity()
  for i in range(3):
    # combine the rotational part by matrix multiplication
    for j in range(3):
      val = 0.0
      for k in range(3):
         val += matrix_elem(a, k, i) * matrix_elem(b, j, k)
      matrix_elem(c, j, i, val)
    val = matrix_elem(a, 3, i)
    # combine the translational part by vector rotation
    for k in range(3):
      val += matrix_elem(a, k, i) * matrix_elem(b, 3, k)
    matrix_elem(c, 3, i, val)
  return c



def is_similar_mag(a, b, small=1E-4):
  """
  Evaluates similar magnitudes to within small.
  """
  return abs(a-b) <= small


def is_similar_vector(a, b, small=1E-4):
  """
  Evaluates similar matrixes through matrix components.
  """
  for a_x, b_x in zip(a, b):
    if not is_similar_mag(a_x, b_x, small):
      return False
  return True


def is_similar_matrix(a, b, small=1E-4):
  """
  Evaluates similar matrixes through matrix components.
  """
  for i in range(0, 3):
    for j in range(0, 3):
      if not is_similar_mag(matrix_elem(a,i,j), matrix_elem(b,i,j), small):
        return False
  return True


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
  Returns component of v perpendicular to axis.
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
  Returns distance between two points
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




