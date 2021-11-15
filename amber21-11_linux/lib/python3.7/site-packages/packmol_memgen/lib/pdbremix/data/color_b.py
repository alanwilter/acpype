# rlc color_b.py version 8.0
# Copyright (c) 2004 Robert L. Campbell
# added user defined colors for 3 color ramp- Mark A. Wall

import colorsys,sys,re
from pymol import cmd

# main function called from within PyMOL
def color_b(selection='all',item='b',mode='hist',gradient='bgr',nbins=11,sat=1.,value=1.,minimum='',maximum='',user_rgb=''):
  """
    
  AUTHOR 

    Robert L. Campbell with enhancements from James Stroud

  USAGE

    color_b selection='sel',item='b', 'q', 'partial_charge' or 'formal_charge'
      gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb' or
      'rw' or 'wr' or 'gw' or 'wg' or 'bw' or wb' or 'gy' or 'yg' or 
      'gray' or 'reversegray' or 'user'

      mode='hist' or 'ramp' (default is 'hist')

      [minimum=''],[maximum=20.],

      nbins=11, sat=1.0, value=1.0, 
      
      user_rgb = '(r1,g1,b1,r2,g2,b2,r3,g3,b3') [for use with gradient=user]

      The "item" argument allows specifying 'b', 'q', 'partial_charge'
      or 'formal_charge'as the item to color on.  The "color_q" function
      is really just the same as "color_b item=q".

      This function allows coloring of a selection as a function of
      B-value or occupancy, following a gradient of colours.  The
      gradients can be:

      'bgr': blue -> green   -> red
      'rgb': red  -> green   -> blue
      'bwr': blue -> white   -> red
      'rwb': red  -> white   -> blue
      'bmr': blue -> magenta -> red
      'rmb': red  -> magenta -> blue
      'rw' : red -> white
      'wr' : white -> red
      'gw' : green -> white
      'wg' : white -> green
      'bw' : blue -> white
      'wb' : white -> blue
      'gy' : green -> yellow
      'yg' : yellow -> green
      'gray' : black -> white
      'reversegray' : white -> black
      'user' : user defined in this script

      ('rainbow' and 'reverserainbow' can be used as synonyms for 
      'bgr' and 'rgb' respectively and 'grey' can be used as a synonym for 'gray').

      User-defined gradients are entered on the command line in
      parentheses as either integers between 0 and 255 or floats between
      0 and 1.  If any one value is larger than 1, then it is assumed
      that all are being entered as integers between 0 and 255.  Hence one can type:

      color_b selection, gradient=user, user_rgb=(0,0,1, 0,.5,1., 1.,.5,0.)

        or

      color_b selection, gradient=user, user_rgb=(0,0,255, 0,128,255, 255,128,0.)

      The division of B-value ranges can in either of two modes: 'hist' or
      'ramp'. 'hist' is like a histogram (equal-sized B-value increments
      leading to unequal numbers of atoms in each bin). 'ramp' as a ramp
      of B-value ranges with the ranges chosen to provide an equal number
      of atoms in each group.

      You can also specify the lower or upper limits of the data used to determine
      the color bins (minimum,maximum). e.g. color_b my_molecule, minimum=15., maximum=25.

      You can also specify the saturation and value (i.e. the "s" and "v"
      in the "HSV" color scheme) to be used for the gradient. The defaults
      are 1.0 for both "sat" and "value".

      In the case of the gray scale gradients, "sat" sets the minimum intensity 
      (normally black) and "value" sets the maximum (normally white)

    usage:
      from within PyMOL do "run color_b.py" to load the function definition.  
      Then you can use for example:

          color_b (c. a | c. b),mode=ramp,gradient=bwr,nbins=30,sat=.5, value=1.

      to color chains A and B with the Blue-White-Red gradient in 30 colors of equal 
      numbers of atoms in each color.
  """

  nbins=int(nbins)
  sat=float(sat)
  value=float(value)
# make sure sat and value are in the range 0-1.0
  sat = min(sat, 1.0)
  sat = max(sat, 0.0)
  value = min(value, 1.0)
  value = max(value, 0.0)
  if gradient == 'user' and user_rgb == '':
    user_rgb = '50,50,195, 245,245,20, 255,20,20'

# make sure lowercase
  gradient.lower()
  mode.lower()

# Sanity checking
  if nbins == 1:
    print("\n     WARNING: You specified nbins=1, which doesn't make sense...resetting nbins=11\n")
    nbins=11

  if mode not in ('hist','ramp'):
    print("\n     WARNING: Unknown mode ",mode, "    ----->   Nothing done.\n")
    return
  elif gradient not in ('bgr','rgb','rainbow','reverserainbow','bwr','rwb','user',
                        'bmr','rmb','rw','wr','gw','wg','bw','wb','gy','yg','gray','grey','reversegray','reversegrey','user'):
    print("\n     WARNING: Unknown gradient: ",gradient, "    ----->   Nothing done.\n")
    return

  print("MODE, GRADIENT, NBINS:", mode,gradient, nbins)

# get list of B-factors from selection
  m = cmd.get_model(selection)
  sel = []
  b_list = []

  if len(m.atom) == 0:
    print("Sorry, no atoms selected")

  else:
    if item == 'b':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].b)
    elif item == 'q':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].q)

    elif item == 'partial_charge':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].partial_charge)

    elif item == 'formal_charge':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].formal_charge)

    else:
      print("Not configured to work on item %s" % item)
      return

    max_b = max(b_list)
    min_b = min(b_list)
    print("Minimum and Maximum B-values: ", min_b, max_b)

    if mode == 'ramp':
      # color in bins of equal numbers of atoms
      b_list.sort()

      # subtract 0.1 from the lowest B in order to ensure that the single
      # atom with the lowest B value doesn't get omitted
#      b_list[0] = b_list[0] - 0.1

      bin_num = int(len(b_list)/nbins)
#      sel.append(selection + " and (b < " + str(b_list[bin_num]) + " or b = " + str(b_list[bin_num]) + ")")
      sel.append(selection + " and (%s < %4.4g" % (item,b_list[bin_num]) + " or %s = %4.4g" % (item,b_list[bin_num]) + ")")
      for j in range(1,nbins):
#        sel.append(selection + " and b > " + str(b_list[j*bin_num]))
        sel.append(selection + " and %s > %4.4g" % (item,b_list[j*bin_num]))
        #print "Color select: ",sel[j]

    elif mode == 'hist':

# check if minimum or maximum was specified and use the entered values
      if minimum != '':
        min_b = float(minimum)
      if maximum != '':
        max_b = float(maximum)
      # histogram:
      # color in bins of equal B-value ranges
      # subtract 0.1 from the lowest B in order to ensure that the single
      # atom with the lowest B value doesn't get omitted
      bin_width = (max_b - min_b)/nbins
      sel.append(selection + " and (%s < %4.4g" % (item,min_b + bin_width) + " or %s = %4.4g" % (item,min_b + bin_width) + ")")
      for j in range(1,nbins):
        sel.append(selection + " and %s > %4.4g" % (item,min_b + j*bin_width))
        #print "Color select: ",sel[j]

# call the function to create the gradient which returns a list of colours
    colours = make_gradient(sel,gradient,nbins,sat,value, user_rgb)

# do the colouring now
    for j in range(nbins):
      print("Color select: ",sel[j])
      cmd.color(colours[j],sel[j])

def color_q(selection="all",mode="hist",gradient="bgr",nbins=11,sat=1.,value=1.,minimum='',maximum=''):
  """
    
  USAGE

    color_q(selection,gradient,mode,nbins,sat,value) ='sel',
      gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb' 
      'rw' or 'wr','gw' or 'wg' or 'bw' or 'wb' or 'gy' or 'yg' or 'gray' or 'reversegray' or 'user'
      mode='hist' or 'ramp', q0=0.,q1=1.0,
      nbins=11, sat=1.0, value=1.0)

      This function allows coloring of a selection as a function of
      occupancy.  See color_b for details.
  """
  item='q'
  color_b(selection,item,mode,gradient,nbins,sat,value,minimum,maximum)

# function for creating the gradient
def make_gradient(sel,gradient,nbins,sat,value,user_rgb):
  if gradient == 'bgr' or gradient == 'rainbow':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # must append the str(sel[j]) to the color name so that it is unique 
      # for the selection
      coldesc.append('col' + str(j))
      # coldesc.append('col' + str(sel[j]) + str(j))

      # create colors using hsv scale (fractional) starting at blue(.6666667) 
      # through red(0.00000) in intervals of .6666667/(nbins -1) (the "nbins-1" 
      # ensures that the last color is, in fact, red (0)
      # rewrote this to use the colorsys module to convert hsv to rgb
      hsv = (colorsys.TWO_THIRD - colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
      #convert to rgb and append to color list
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])
      # print "defined as ", str(sel[j])

  elif gradient == 'user':
    # --------------------------------------------
    #  Modified color ramp by Mark Wall 2007.07.20
    # --------------------------------------------
    # assign 3-color ramp values (rgb 0-255)
    # could easily assign RGB color values between 0.0 and 1.0 below
    # !!! Black must be at least 1,0,0 or div by zero error !!!
    #
    #r1, g1, b1 = 255, 255, 225   # low   white
    #r1, g1, b1 = 170, 170, 170   # low   gray
    user_rgb = re.compile('[\[\](){}]').sub('',user_rgb)
    user_rgb_fract = 0
    try:
      r1,g1,b1, r2,g2,b2, r3,g3,b3 = list(map(int,user_rgb.split(',')))
    except ValueError:
      r1,g1,b1, r2,g2,b2, r3,g3,b3 = list(map(float,user_rgb.split(',')))
      user_rgb_fract = 1
      
    print('user_rgb', r1,g1,b1, r2,g2,b2, r3,g3,b3)

#    r1, g1, b1 = 50, 50, 195     # low   med blue
#    r2, g2, b2 = 245, 245, 20    # mid   yellow
#    r3, g3, b3 = 255, 20, 20     # high  red
    #
    #r1, g1, b1 = 255, 20, 20     # low   red
    #r2, g2, b2 = 150, 150, 20    # mid   yellow
    #r3, g3, b3 = 20, 20, 195     # high   med blue
    #
    #r1, g1, b1 = 1, 0, 0         # low   black
    #r2, g2, b2 = 155, 155, 155   # mid  gray
    #r3, g3, b3 = 255, 255, 255   # high  white
    #
    #r1, g1, b1 = 0, 50, 200      # low  blue
    #r2, g2, b2 = 1, 0, 0         # mid   black
    #r3, g3, b3 = 255, 255, 20    # high  yellow
    #
    #r1, g1, b1 = 0, 0, 1         # low   black
    #r2, g2, b2 = 200, 0, 0       # mid  red
    #r3, g3, b3 = 255, 255, 0     # high  yellow
    #
    #r1, g1, b1 = 180, 170, 170   # low   gray
    #r2, g2, b2 = 250, 90, 40     # mid  orange
    #r3, g3, b3 = 255, 255, 0     # high  yellow
    #
    #r1, g1, b1 = 235, 255, 255   # low   white
    #r2, g2, b2 = 55, 255, 255    # mid   cyan
    #r3, g3, b3 = 0, 0, 180       # high  blue
    #
    # change color values to fractions
    # 
    if max(r1,g1,b1,r2,g2,b2,r3,g3,b3) > 1:
      r1, g1, b1 = float(r1)/255.0, float(g1)/255.0, float(b1)/255.0
      r2, g2, b2 = float(r2)/255.0, float(g2)/255.0, float(b2)/255.0 
      r3, g3, b3 = float(r3)/255.0, float(g3)/255.0, float(b3)/255.0

    col=[]
    coldesc=[]
#    print "r1,g1,b1, r2,g2,b2, r3,g3,b3", r1,g1,b1, r2,g2,b2, r3,g3,b3
    for j in range(nbins/2):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from low to mid 
      rgb = [r1*((float(nbins)-float(j)*2.0)/float(nbins))+r2*(float(j)*2.0/float(nbins)), \
             g1*((float(nbins)-float(j)*2.0)/float(nbins))+g2*(float(j)*2.0/float(nbins)), \
             b1*((float(nbins)-float(j)*2.0)/float(nbins))+b2*(float(j)*2.0/float(nbins))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])
      print(j,"rgb: %4.3f %4.3f %4.3f"% rgb)

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

    for j in range(nbins/2,nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from mid to high
      rgb = [r2*((float(nbins)-((float(j+1)-float(nbins)/2.0)*2.0))/float(nbins))+r3*(((float(j+1)-float(nbins)/2.0)*2.0)/float(nbins)), \
             g2*((float(nbins)-((float(j+1)-float(nbins)/2.0)*2.0))/float(nbins))+g3*(((float(j+1)-float(nbins)/2.0)*2.0)/float(nbins)), \
             b2*((float(nbins)-((float(j+1)-float(nbins)/2.0)*2.0))/float(nbins))+b3*(((float(j+1)-float(nbins)/2.0)*2.0)/float(nbins))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])
      print(j,"rgb: %4.3f %4.3f %4.3f"% rgb)

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'rgb' or gradient == 'reverserainbow':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # must append the str(sel[j]) to the color name so that it is unique 
      # for the selection
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))

      # create colors using hsv scale (fractional) starting at red(.00000) 
      # through blue(0.66667) in intervals of .6666667/(nbins -1) (the "nbins-1" 
      # ensures that the last color is, in fact, red (0)
      # rewrote this to use the colorsys module to convert hsv to rgb
      hsv = (colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
      #convert to rgb and append to color list
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'bmr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from blue through magenta to red
      rgb = [min(1.0, float(j)*2/(nbins-1)), 0.0, min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'rmb':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from red through magenta to blue
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), 0.0, min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'rw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from red through white
      rgb = [1.0, float(j)/(nbins-1), float(j)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'wr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from white through red 
      rgb = [1.0, float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'gw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from green through white
      rgb = [float(j)/(nbins-1), 1.0, float(j)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'wg':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from white through green 
      rgb = [float(nbins-j-1)/(nbins-1), 1.0, float(nbins-j-1)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'bw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from blue through white
      rgb = [float(j)/(nbins-1), float(j)/(nbins-1), 1.0 ]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'wb':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from blue through white
      rgb = [float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1), 1.0 ]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'wr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from white through blue 
      rgb = [float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1), 1.0]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'gy':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from green through yellow
      rgb = [float(j)/(nbins-1), 1.0, 0.]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'yg':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from green through yellow
      rgb = [float(nbins-j-1)/(nbins-1), 1.0, 0.]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'bwr':
    col=[]
    coldesc=[]
    for j in range(nbins/2):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from blue to white 
      rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

    for j in range(nbins/2,nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from white to red
      rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'rwb':
    col=[]
    coldesc=[]
    for j in range(nbins/2):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from red to white
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

    for j in range(nbins/2,nbins):
      coldesc.append('col' + str(j))
      # coldesc.append('col' + str(sel[j]) + str(j))
      # create colors in a gradient from white to blue
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'gray' or gradient == 'grey':
# if it is "gray" then sat must be 0!
    sat = 0.0
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient of grays from "sat" to "value"

      hsv = [0, 0, sat + (value-sat)*float(j)/(nbins-1)]
#      hsv[1] = hsv[1]*sat
#      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'reversegray' or gradient == 'reversegrey':
# if it is "gray" then sat must be 0!
    sat = 0.0
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient of grays from "sat" to "value"

      hsv = [0, 0, value - (value-sat)*float(j)/(nbins-1)]
#      hsv[1] = hsv[1]*sat
#      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])


# return the gradient as a list of colors named by their index (i.e. col0,col1,col2,col3,...)
  return coldesc

# allow calling without parentheses: color_hist_b [selection=], [mode= ],[gradient= ],[nbins= ]
cmd.extend("color_b",color_b)
cmd.extend("color_q",color_q)

