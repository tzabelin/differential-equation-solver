set terminal pngcairo enhanced font 'Verdana,10'
set output '3dplot_multiple_trajectories_Y.png'

# Set the labels
set xlabel 'X-axis'
set ylabel 'Y-axis'
set zlabel 'Z-axis'

# Set the title
set title '3D Plot of 55 Trajectories'

# Enable grid
set grid

# Set 3D view
set view 90,0,1,1
splot \
  'output.dat' using 1:2:3 with lines title '', \
  'output.dat' using 4:5:6 with lines title '', \
  'output.dat' using 7:8:9 with lines title '', \
  'output.dat' using 10:11:12 with lines title '', \
  'output.dat' using 13:14:15 with lines title '', \
  'output.dat' using 16:17:18 with lines title '', \
  'output.dat' using 19:20:21 with lines title '', \
  'output.dat' using 22:23:24 with lines title '', \
  'output.dat' using 25:26:27 with lines title '', \
  'output.dat' using 28:29:30 with lines title '', \
  'output.dat' using 31:32:33 with lines title '', \
  'output.dat' using 34:35:36 with lines title '', \
  'output.dat' using 37:38:39 with lines title '', \
  'output.dat' using 40:41:42 with lines title '', \
  'output.dat' using 43:44:45 with lines title '', \
  'output.dat' using 46:47:48 with lines title '', \
  'output.dat' using 49:50:51 with lines title '', \
  'output.dat' using 52:53:54 with lines title '', \
  'output.dat' using 55:56:57 with lines title '', \
  'output.dat' using 58:59:60 with lines title '', \
  'output.dat' using 61:62:63 with lines title '', \
  'output.dat' using 64:65:66 with lines title '', \
  'output.dat' using 67:68:69 with lines title '', \
  'output.dat' using 70:71:72 with lines title '', \
  'output.dat' using 73:74:75 with lines title '', \
  'output.dat' using 76:77:78 with lines title '', \
  'output.dat' using 79:80:81 with lines title '', \
  'output.dat' using 82:83:84 with lines title '', \
  'output.dat' using 85:86:87 with lines title '', \
  'output.dat' using 88:89:90 with lines title '', \
  'output.dat' using 91:92:93 with lines title '', \
  'output.dat' using 94:95:96 with lines title '', \
  'output.dat' using 97:98:99 with lines title '', \
  'output.dat' using 100:101:102 with lines title '', \
  'output.dat' using 103:104:105 with lines title '', \
  'output.dat' using 106:107:108 with lines title '', \
  'output.dat' using 109:110:111 with lines title '', \
  'output.dat' using 112:113:114 with lines title '', \
  'output.dat' using 115:116:117 with lines title '', \
  'output.dat' using 118:119:120 with lines title '', \
  'output.dat' using 121:122:123 with lines title '', \
  'output.dat' using 124:125:126 with lines title '', \
  'output.dat' using 127:128:129 with lines title '', \
  'output.dat' using 130:131:132 with lines title '', \
  'output.dat' using 133:134:135 with lines title '', \
  'output.dat' using 136:137:138 with lines title '', \
  'output.dat' using 139:140:141 with lines title '', \
  'output.dat' using 142:143:144 with lines title '', \
  'output.dat' using 145:146:147 with lines title '', \
  'output.dat' using 148:149:150 with lines title '', \
  'output.dat' using 151:152:153 with lines title '', \
  'output.dat' using 154:155:156 with lines title '', \
  'output.dat' using 157:158:159 with lines title '', \
  'output.dat' using 160:161:162 with lines title '', \
  'output.dat' using 163:164:165 with lines title ''

