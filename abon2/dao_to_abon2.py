daolist=['67494',
	 '89531',
	 '101395',
	 '234779',
	 '234816',
	 '234822',
	 '234932',
	 '244523',
	 '244652',
	 '244694',
	 '244819',
	 '244853',
	 '256289',
	 '256322',
	 '392942',
	 '402322',
	 '402370',
	 '402384',
	 '402386',
	 '412759',
	 '423375',
	 '554709'
	]
for i in range(len(daolist)):
	filnam=open('H11/'+daolist[i]+'.daospec', 'r')
	lbd=[]
	spc=[]
	ew=[]
	chiex=[]
	loggf=[]
	ref=[]
	c6=[]
	err=[]
	qf=[]
	for cols in (raw.strip().split() for raw in filnam):
		if (cols[0][0] != 'R') and (len(cols) > 5):
			lbd.append(float(cols[5]))
			spc.append(cols[6])
			ew.append(float(cols[2]))
			chiex.append(float(cols[7]))
			ref.append(cols[8])
			loggf.append(float(cols[9]))
			c6.append(cols[10])
			err.append(float(cols[3]))
			qf.append(float(cols[4]))
	filnam.close()
	filnam=open('H12/'+daolist[i]+'.daospec', 'r')
	for cols in (raw.strip().split() for raw in filnam):
		if (cols[0][0] != 'R') and (len(cols) > 5) and (cols[6] != 'FE2'):
			lbd.append(float(cols[5]))
			spc.append(cols[6])
			ew.append(float(cols[2]))
			chiex.append(float(cols[7]))
			ref.append(cols[8])
			loggf.append(float(cols[9]))
			c6.append(cols[10])
			err.append(float(cols[3]))
			qf.append(float(cols[4]))
	filnam.close()
	outfil=open('raies_'+daolist[i]+'.dat', 'w')
	outfil.write('####\n')
	outfil.write('FORMAT(F8.3,1X,A2,I1,F9.6,I4,F7.3,2X,e9.3,3X,F6.1)\n')
	outfil.write('lambda   el kiex ref log gf      C6        EW    err     QF\n')
	for j in range(len(lbd)):
		outfil.write('%7.3f%4s%9.6f%4s%7.3f%11s%9.1f%6.1f%8.3f\n' %(lbd[j], spc[j], chiex[j], ref[j], loggf[j], c6[j], ew[j], err[j], qf[j]))
	outfil.close()
