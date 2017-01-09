daolist=['6159',
	 '18025',
	 '18731',
	 '28233',
	 '28326',
	 '28734',
	 '34152'
	]
for i in range(len(daolist)):
	filnam=open(daolist[i]+'_L.daospec', 'r')
	lbd=[]
	ew=[]
	err=[]
	for cols in (raw.strip().split() for raw in filnam):
		if (cols[0][0] != 'R') and (len(cols) > 5):
			lbd.append(cols[5])
			ew.append(cols[2])
			err.append(cols[3])
	filnam.close()
	filnam=open(daolist[i]+'_U.daospec', 'r')
	for cols in (raw.strip().split() for raw in filnam):
		if (cols[0][0] != 'R') and (len(cols) > 5):
			lbd.append(cols[5])
			ew.append(cols[2])
			err.append(cols[3])
	filnam.close()
	template=open('model.list', 'r')
	outfil=open('syn_'+daolist[i]+'.list', 'w')
	for cols in (raw.strip().split() for raw in template):
		if (len(cols) == 4):
			outfil.write(cols[0]+cols[1]+' '+cols[2]+' '+cols[3]+'\n')
		elif (len(cols) == 2):
			outfil.write(cols[0]+' '+cols[1]+'\n')
		else:
			flag=False
			for i in range(len(lbd)):
				if cols[0] == lbd[i]:
					outfil.write('  '+cols[0]+'  '+cols[1]+'%8.3f'%(float(cols[2]))+'  '+cols[3]+'%7.1f'%(float(cols[4]))+'  '+cols[5]+' '+cols[6]+' '+cols[7]+'%8.2f'%(float(ew[i]))+'    '+err[i]+'\n')
					flag=True
			if flag == False:
				outfil.write('  '+cols[0]+'  '+cols[1]+'%8.3f'%(float(cols[2]))+'  '+cols[3]+'%7.1f'%(float(cols[4]))+'  '+cols[5]+' '+cols[6]+' '+cols[7]+'%8.2f'%(0.0)+'    '+'0.0'+'\n')
	outfil.close()
	template.close()
