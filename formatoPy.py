from pprint import pprint
from sl3offsets import sl3_ofssets

## event data
jaime_dump = """
Wh:-2 Se:7 St:4
37 79645 0 37 79645 0
37 79414 1 37 79414 1
37 79621 0 37 79621 0
37 79407 1 37 79407 1
-1 -1 -1 35 79410 0
-1 -1 -1 35 79688 1
-1 -1 -1 35 79391 0
-1 -1 -1 -1 -1 -1
PosEmul 155.844 PosFW 156.55
TimeEmul 79326 TimeFW 79338
TanPsiEmul -0.0161133 TanPsiFW -0.0600586
chi2Emul 0.00664062 chi2FW 0.0247168
"""

jaime_dump = """
Wh:0 Se:6 St:3
-1 -1 -1 -1 -1 -1
-1 -1 -1 -1 -1 -1
-1 -1 -1 -1 -1 -1
-1 -1 -1 -1 -1 -1
15 34207 0 14 33939 1
15 33939 1 15 33939 0
14 34276 1 14 33959 1
15 33952 0 15 33952 0
PosEmul 63.075 PosFW 61.9562
TimeEmul 33898 TimeFW 33754
TanPsiEmul 0.185059 TanPsiFW -0.00756836
chi2Emul 0.000693359 chi2FW 0.00769531
"""


data = {}

[wh, se, st] =  [ int(jaime_dump.split('\n')[1].replace(':',' ').split(' ')[i]) for i in [1,3,5]]
data['wh'] = wh
data['se'] = se
data['st'] = st

sl_phi3_x_offset = sl3_ofssets[(wh,se,st)]


hitlines = jaime_dump.split('\n')[2:10]
for (sl, slstart) in [('phi1', 0), ('phi2',4)]:
  data[sl]={}
  data[sl]['valid'] = [ '-' not in hitlines[i]          for i in range(slstart,slstart+4) ]
  data[sl]['wires'] = [ int(hitlines[i].split(' ')[0])  for i in range(slstart,slstart+4) ]
  data[sl]['t0s']   = [ int(hitlines[i].split(' ')[1])  for i in range(slstart,slstart+4) ]
  data[sl]['lat']   = [ int(hitlines[i].split(' ')[2])  for i in range(slstart,slstart+4) ]

for (info,t) in [('Pos',float),('Time',int),('TanPsi',float),('chi2',float)]:
  data[info]={
    'FW'  : t(jaime_dump.split(info+'FW '  )[1].split('\n')[0]),
    'Emul': t(jaime_dump.split(info+'Emul ')[1].split(' ' )[0]),
    }

#t0 = []
#cells = []

sls = ['phi1','phi2']

for sl in range(0,2):
  t0s = data[sls[sl]]['t0s']
  wires = data[sls[sl]]['wires']
  for la in range(0,4):
    if not t0s[la] == -1 : print wh, se, st, sl_phi3_x_offset, 2*sl, la, wires[la], t0s[la]  



