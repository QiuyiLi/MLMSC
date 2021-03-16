import os

names = ['N50000D5.7L4.9C9P1', 'N50000D5.7L4.9C9P2', 'N50000D5L5C9P1', 'N50000D5L5C9P2', 'N50000D4L4C9P1', 'N50000D4L4C9P2', 'N50000D3L3C9P1', 'N50000D3L3C9P2', 'N50000D2L2C9P1', 'N50000D2L2C9P2', 'N50000D1L1C9P1', 'N50000D1L1C9P2', 'N50000D0L0C9P1', 'N50000D0L0C9P2']
f = open('template.slurm','r')
string = f.read()
for name in names:
	fileName = name + '.slurm'
	g = open(fileName,'w')
	g.write(string.replace('N30D2L2',name))
	g.close()
f.close()