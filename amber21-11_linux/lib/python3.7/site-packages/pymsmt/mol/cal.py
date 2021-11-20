"""
This module is used for the basic calculation of a molecule system.
"""
import math
import numpy

#------------------------------------------------------------------------------

def calc_bond(at1, at2):
    bond = math.sqrt((at1[0]-at2[0])**2+(at1[1]- at2[1])**2+(at1[2]-at2[2])**2)
    return bond

##Use cosine law to calculate the angle
def calc_angle(at1, at2, at3):
    d12 = calc_bond(at1, at2)
    d23 = calc_bond(at2, at3)
    d13 = calc_bond(at1, at3)
    tempval = (d23**2+d12**2-d13**2)/(2*d12*d23)
    #print "%32.16f" %tempval
    if (tempval <= -1.0):
        angle = 180.0
    elif (tempval >= 1.0):
        angle = 0.0
    else:
        angle = 180.0*math.acos(tempval)/math.pi
    return angle

##The functions for vector calculation use
def calc_vec_cross(a, b):
    #(a0, a1, a2) * ( b0, b1, b2) = (a1b2 - a2b1, a2b0 - a0b2, a0b1 - a1b0)
    result = ( a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
               a[0] * b[1] - a[1] * b[0])
    return result

def calc_vec_dot(a, b):
    #(x1, y1, z1).(x2, y2, z2) = x1x2 + y1y2 + z1z2
    result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    return result

def calc_vec_min(a, b):
    result = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    return result

def calc_vec_value(a):
    result = math.sqrt((a[0])**2+(a[1])**2+(a[2])**2)
    return result

##Calculate the dihedrals
def calc_dih(at1, at2, at3, at4):

    b1 = calc_vec_min(at2, at1)
    b2 = calc_vec_min(at3, at2)
    b3 = calc_vec_min(at4, at3)

    b12 = calc_vec_cross(b1, b2)
    b23 = calc_vec_cross(b2, b3)
    b2vv = calc_vec_value(b2)

    term1in1 = calc_vec_cross(b12, b23)
    term2in1 = (b2[0]/ b2vv, b2[1]/b2vv, b2[2]/b2vv)

    term2 = calc_vec_dot(b12, b23)

    term1 = calc_vec_dot(term1in1, term2in1)

    dih = math.atan2(term1, term2)
    dih = 180 * dih / math.pi

    return dih

def calc_dih_fix(bond12, bond23, bond34, ang123, ang234):
    nar1 = bond12 * math.sin(math.radians(180.0-ang123))
    nar2 = bond34 * math.sin(math.radians(180.0-ang234))
    nar3 = bond34 * math.cos(math.radians(180.0-ang234))
    nar4 = bond12 * math.cos(math.radians(180.0-ang123))
    tmp = (nar1**2 + nar2**2 - dis**2 + (nar3 + bond23 + nar4)**2)/\
          (2 * nar1 * nar2)
    if (tmp >= -1.0) and (tmp <= 1.0):
        dih = math.acos(tmp)
        return dih
    else:
        return "Error"

def get_angles(metcrd, crds):

    angles = []
    for i in range(0, len(crds)):
        crdi = crds[i]
        for j in range(i+1, len(crds)):
            crdj = crds[j]
            angle = calc_angle(crdi, metcrd, crdj)
            angles.append(angle)

    angles.sort()
    return angles

def det_geo(crds):

    metcrd = crds[0]
    crds = crds[1::2]
    angles = get_angles(metcrd, crds)

    if len(crds) == 1: #1-Coordinated, no angle
        angrms = 0.0
        angrms = round(angrms, 3)
        return '1', angrms
    elif len(crds) == 2:  #2-Coordinated, 1 angle
        #Linear
        angrms = [i-180.0 for i in angles]
        angrms = [abs(i) for i in angrms]
        angrms = numpy.sqrt(numpy.average(angrms))
        angrms = round(angrms, 3)
        return '2Ln', angrms
    elif len(crds) == 3: #3-Coordinated, 3 angles
        #Triangle
        angrms = [i-120.0 for i in angles]
        angrms = [abs(i) for i in angrms]
        angrms = numpy.sqrt(numpy.average(angrms))
        angrms = round(angrms, 3)
        return '3Tr', angrms
    elif len(crds) == 4: #4-Coordinated, 6 angles
        #1. Tetrahedral
        angrms1 = [i-109.5 for i in angles]
        angrms1 = [abs(i) for i in angrms1]
        angrms1 = numpy.sqrt(numpy.average(angrms1))
        #2. Square Planar
        angrms21 = [i-90.0 for i in angles[0:4]]
        angrms22 = [i-180.0 for i in angles[4:6]]
        angrms2 =  angrms21 + angrms22
        angrms2 = [abs(i) for i in angrms2]
        angrms2 = numpy.sqrt(numpy.average(angrms2))
        if (angrms1 < angrms2):
            return '4Te', angrms1
        else:
            return '4Sq', angrms2
    elif len(crds) == 5: #5-Coordinated, 10 angles
        angrmss = []
        #1. Trigonal Bipyramid
        angrms11 = [i-90.0 for i in angles[0:6]]
        angrms12 = [i-120.0 for i in angles[6:9]]
        angrms13 = [i-180.0 for i in angles[9:10]]
        angrms1 = angrms11 + angrms12 + angrms13
        angrms1 = [abs(i) for i in angrms1]
        angrms1 = numpy.sqrt(numpy.average(angrms1))
        angrmss.append(angrms1)

        #2. Square Pyramid
        tmpangs = []
        angrms2s = []
        #try each opptunity
        for i in range(0, len(crds)):
            topcrd = crds[i]
            othcrd = crds[:i] + crds[i+1:]
            for j in range(0, len(othcrd)):
                angle = calc_angle(topcrd, metcrd, othcrd[j])
                tmpangs.append(angle)
            avgang1 = numpy.average(tmpangs)
            avgang2 = 360.0 - 2 * avgang1
            avgang3 = 2 * math.asin(1.0/math.sqrt(2.0) * (math.sin(180.0 - avgang1)))
            avgangs = []
            for k in range(0, 4):
                avgangs.append(avgang1)
                avgangs.append(avgang3)
            for k in range(0, 2):
                avgangs.append(avgang2)
            avgangs.sort()
            angrms2 = [abs(angles[i]-avgangs[i]) for i in range(0, 10)]
            angrms2 = numpy.sqrt(numpy.average(angrms2))
            angrms2s.append(angrms2)
        angrms2 = min(angrms2s)
        angrmss.append(angrms2)

        #3. Tetrahedral with Nonbonded interaction
        angrms3s = []
        #Try each opportunity
        for i in range(0, len(crds)):
            tmpcrds = crds[0:i] + crds[i+1:]
            tmpangs = get_angles(metcrd, tmpcrds)
            angrms3 = [i-109.5 for i in tmpangs]
            angrms3 = [abs(i) for i in angrms3]
            angrms3 = numpy.sqrt(numpy.average(angrms3))
            angrms3s.append(angrms3)
        angrms3 = min(angrms3s)
        angrmss.append(angrms3)

        if min(angrmss) == angrms1:
            angrms1 = round(angrms1, 3)
            return '5Tp', angrms1
        elif min(angrmss) == angrms2:
            angrms2 = round(angrms2, 3)
            return '5Sp', angrms2
        elif min(angrmss) == angrms3:
            angrms3 = round(angrms3, 3)
            return '5Tn', angrms3
    elif len(crds) == 6: #6-Coordinated, 15 angles
        #Octohedral
        angrms11 = [i-90.0 for i in angles[0:12]]
        angrms12 = [i-180.0 for i in angles[12:15]]
        angrms = angrms11 + angrms12
        angrms = [abs(i) for i in angrms]
        angrms = numpy.sqrt(numpy.average(angrms))
        angrms = round(angrms, 3)
        return '6Oc', angrms
    elif len(crds) == 7: #7-Coordinated, 21 angles
        angrms11 = [i-72.0 for i in angles[0:5]]
        angrms12 = [i-90.0 for i in angles[5:14]]
        angrms13 = [i-144.0 for i in angles[14:20]]
        angrms14 = angles[20]-180.0
        angrms = angrms11 + angrms12 + angrms13
        angrms.append(angrms14)
        angrms = [abs(i) for i in angrms]
        angrms = numpy.sqrt(numpy.average(angrms))
        angrms = round(angrms, 3)
        return '7Bt', angrms
    elif len(crds) == 8: #8-Coordinated, 28 angles
        angrms11 = [i-90.0 for i in angles[0:24]]
        angrms12 = [i-180.0 for i in angles[24:28]]
        angrms = angrms11 + angrms12
        angrms = [abs(i) for i in angrms]
        angrms = numpy.sqrt(numpy.average(angrms))
        angrms = round(angrms, 3)
        return '8Bt', angrms
