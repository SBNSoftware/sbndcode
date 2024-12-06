#!/usr/bin/env python

rayzx_list = []
rayzx_obj_list = []

def print_ray(head, tail, halfD=1.8):
  print '''  {
    head: wc.point(%.3f, %.3f, %.3f, wc.m),
    tail: wc.point(%.3f, %.3f, %.3f, wc.m),  
  },
  ''' %(head[0], head[1], head[2], tail[0], tail[1], tail[2])

  print '''  {
    head: wc.point(%.3f, %.3f, %.3f, wc.m),
    tail: wc.point(%.3f, %.3f, %.3f, wc.m),  
  },
  ''' %(head[0] + halfD*2, head[1], head[2], tail[0] + halfD*2, tail[1], tail[2])
  rayzx_list.append((head[2], head[0], tail[2], tail[0]))
  rayzx_list.append((head[2], head[0] + halfD*2, tail[2], tail[0]+ halfD*2))


import sys, math
argc = len(sys.argv)
if argc!=2:
  print "Usage example: python gen.py 45 # unit in degree" 

deg = math.radians(1)
theta = float(sys.argv[1]) *deg

halfZ = 3.45 # meter
halfD = 1.8 # half of drift length
halfY = 3.0


print '''local wc = import 'wirecell.jsonnet';
{
  rays: [
'''


head = None
tail = None
if halfZ * math.tan(theta) < halfD:
  head = [-halfD - halfZ*math.tan(theta), halfY, 0.0]
  tail = [-halfD + halfZ*math.tan(theta), halfY, halfZ*2]
  print_ray(head,tail)

else:
  head = [-halfD*2, halfY, halfZ - halfD/math.tan(theta)]
  tail = [     0.0, halfY, halfZ + halfD/math.tan(theta)]

  half_ray_span = halfD/math.tan(theta)
  half_extra_nrays = int( (halfZ-half_ray_span) / (2*half_ray_span + 2*0.05) ) + 1 # 0.05 gaps between tracks
  if half_extra_nrays<1: half_extra_nrays = 1

  print_ray(head,tail)
  for i in range(half_extra_nrays):
    _head = head[:]  
    _tail = tail[:]  
    shiftZ = (i+1) * (half_ray_span*2 + 0.05) 
    _head[2] -= shiftZ
    _tail[2] -= shiftZ
    # print "shiftZ: ", shiftZ
    # print "_head[0]: ", _head[0]
    if _head[2] < 0.0: 
      _head[0] -= math.tan(theta) * _head[2]
      _head[2] = 0.0
  
    head_ = head[:]  
    tail_ = tail[:]  
    head_[2] += shiftZ
    tail_[2] += shiftZ
    if tail_[2] > (halfZ*2):
      tail_[0] -= math.tan(theta) * (tail_[2] - halfZ*2)
      tail_[2] = halfZ*2

    print_ray(_head, _tail)
    print_ray(head_, tail_)


print '''  ],
}
'''


from ROOT import *
c1 = TCanvas("c1","c1",800,600)
h = TH2F("h","", 100,0-0.5,6.9+0.5,100,-3.6-0.5,3.6+0.5)
h.Draw()

# rayzx_list = [(1.65, -3.6, 5.25, 0.0), (1.65, 0.0, 5.25, 3.6)]
for zx in rayzx_list:
  ray1 = TLine(zx[0], zx[1], zx[2], zx[3])
  rayzx_obj_list.append(ray1)
  ray1.Draw()
