import math
import numpy as np
import matplotlib.pyplot as plt

def mmq_adj(x,y):
	n=len(x)
	A=[[n,0],[0,0]]
	B=[0,0]
	for i in range(n):
		A[0][1]=A[0][1]+x[i]
		A[1][0]=A[1][0]+x[i]
		A[1][1]=A[1][1]+(x[i]**2)
		B[0]=B[0]+y[i]
		B[1]=B[1]+x[i]*y[i]
	p=(A[1][0]/A[0][0])
	A[1][1]=A[1][1]-(p*A[0][1])
	B[1]=B[1]-(p*B[0])
	slope=B[1]/A[1][1]
	inter=(B[0]-A[0][1]*slope)/A[0][0]
	return slope,inter

def colorplot(X, Y, lbd):
	for n, i in enumerate(X):
		if lbd[n] < 4400:
			plt.plot(X[n], Y[n], color='#800080', marker='.')
		elif lbd[n] > 4400 and lbd[n] < 4850:
			plt.plot(X[n], Y[n], color='#0000FF', marker='.')
		elif lbd[n] > 4850 and lbd[n] < 5000:
			plt.plot(X[n], Y[n], color='#00FFFF', marker='.')
		elif lbd[n] > 5000 and lbd[n] < 5650:
			plt.plot(X[n], Y[n], color='#00DD00', marker='.')
		elif lbd[n] > 5650 and lbd[n] < 5900:
			plt.plot(X[n], Y[n], color='#F6F600', marker='.')
		elif lbd[n] > 5900 and lbd[n] < 6250:
			plt.plot(X[n], Y[n], color='#FF8C00', marker='.')
		elif lbd[n] > 6250 and lbd[n] < 7400:
			plt.plot(X[n], Y[n], color='#FF0000', marker='.')

solar=7.5

arq=raw_input()
mult=input()
with open(arq) as f:
	lines = f.readlines()[2:]
el=[]
ion=[]
chiex=[]
lbd=[]
rdew=[]
abon=[]
for cols in (raw.strip().split() for raw in lines):
	el.append(cols[1][1:3])
	ion.append(cols[1][3:4])
	chiex.append(float(cols[2]))
	rdew.append(math.log10((float(cols[6])/1000)/float(cols[0])))
	lbd.append(float(cols[0]))
	abon.append(float(cols[10])-solar)
nlines=len(el)
elion=[]
elion.append([])
elm_chiex=[]
elm_chiex.append([])
elm_rdew=[]
elm_rdew.append([])
elm_lbd=[]
elm_lbd.append([])
elm_abon=[]
elm_abon.append([])
index=0
elements=[]
elements.append(el[1])
for i in range(nlines):
	if i == (nlines-1):
		elm_chiex[index].append(chiex[i])
		elm_rdew[index].append(rdew[i])
		elm_lbd[index].append(lbd[i])
		elm_abon[index].append(abon[i])
		elion[index].append(ion[i])
	elif el[i] == el[i+1]:
		elm_chiex[index].append(chiex[i])
		elm_rdew[index].append(rdew[i])
		elm_lbd[index].append(lbd[i])
		elm_abon[index].append(abon[i])
		elion[index].append(ion[i])
	else:
		elm_chiex.append([])
		elm_rdew.append([])
		elm_lbd.append([])
		elm_abon.append([])
		elion.append([])
		index += 1
		elements.append(el[i+1])
chiex_i=[]
chiex_ii=[]
rdew_i=[]
lbd_i=[]
rdrdew_ii=[]
lbd_ii=[]
abon_i=[]
abon_ii=[]
nelem=len(elion)
for i in range(nelem):
	chiex_i.append([])
	chiex_ii.append([])
	rdew_i.append([])
	rdrdew_ii.append([])
	lbd_i.append([])
	lbd_ii.append([])
	abon_i.append([])
	abon_ii.append([])
	for j in range(len(elion[i])):
		if elion[i][j] == '1':
			chiex_i[i].append(elm_chiex[i][j])
			rdew_i[i].append(elm_rdew[i][j])
			lbd_i[i].append(elm_lbd[i][j])
			abon_i[i].append(elm_abon[i][j])
		else:
			chiex_ii[i].append(elm_chiex[i][j])
			rdrdew_ii[i].append(elm_rdew[i][j])
			lbd_ii[i].append(elm_lbd[i][j])
			abon_ii[i].append(elm_abon[i][j])
achiex_i=[0 for i in range(nelem)]
achiex_ii=[0 for i in range(nelem)]
bchiex_i=[0 for i in range(nelem)]
bchiex_ii=[0 for i in range(nelem)]
ardew_i=[0 for i in range(nelem)]
ardrdew_ii=[0 for i in range(nelem)]
brdew_i=[0 for i in range(nelem)]
brdrdew_ii=[0 for i in range(nelem)]
abondi=[]
abondii=[]
for n in range(nelem):
	achiex_i[n],bchiex_i[n]=mmq_adj(chiex_i[n],abon_i[n])
	ardew_i[n],brdew_i[n]=mmq_adj(rdew_i[n],abon_i[n])
	abondi.append(sum(abon_i[n])/len(abon_i[n]))
	abondii.append(sum(abon_ii[n])/len(abon_ii[n]))
residuals1=[]
residuals2=[]
rej_abon=[]
rej_chiex=[]
rej_ew=[]
rej_lbd=[]
n_abon=[]
n_chiex=[]
n_ew=[]
n_lbd=[]
for i in range(len(abon_i[0])):
	residuals1.append(math.fabs((chiex_i[0][i]*achiex_i[0]+bchiex_i[0])-(abon_i[0][i])))
	residuals2.append(math.fabs((rdew_i[0][i]*ardew_i[0]+brdew_i[0])-(abon_i[0][i])))

rej=mult*max([np.std(residuals1), np.std(residuals2)])

for i in range(len(residuals1)):
	if residuals1[i] > rej:
		rej_abon.append(abon_i[0][i])
		rej_chiex.append(chiex_i[0][i])
		rej_ew.append(rdew_i[0][i])
		rej_lbd.append(lbd_i[0][i])
	else:
		n_abon.append(abon_i[0][i])
		n_chiex.append(chiex_i[0][i])
		n_ew.append(rdew_i[0][i])
		n_lbd.append(lbd_i[0][i])
###		
achiex_i=[0 for i in range(nelem)]
achiex_ii=[0 for i in range(nelem)]
bchiex_i=[0 for i in range(nelem)]
bchiex_ii=[0 for i in range(nelem)]
ardew_i=[0 for i in range(nelem)]
ardrdew_ii=[0 for i in range(nelem)]
brdew_i=[0 for i in range(nelem)]
brdrdew_ii=[0 for i in range(nelem)]
abondi=[]
for n in range(nelem):
	achiex_i[n],bchiex_i[n]=mmq_adj(n_chiex,n_abon)
	ardew_i[n],brdew_i[n]=mmq_adj(n_ew,n_abon)
	abondi.append(sum(n_abon)/len(n_abon))
###
line1x_i=[-50, 50]
line1y_i=[line1x_i[0]*achiex_i[0]+bchiex_i[0], line1x_i[1]*achiex_i[0]+bchiex_i[0]]
line2x_i=[-10, 200]
line2y_i=[line2x_i[0]*ardew_i[0]+brdew_i[0], line2x_i[1]*ardew_i[0]+brdew_i[0]]
line1x_ii=[-50, 50]
line1y_ii=[abondii[0], abondii[0]]
line2x_ii=[-10, 200]
line2y_ii=[abondii[0], abondii[0]]

with open(arq) as f:
	pars = f.readlines()[1:2]
for cols in (raw.strip().split() for raw in pars):
	stpars='teff='+cols[1]+' logg='+cols[2]+' [Fe/H]='+cols[3]+' vt='+cols[0]

fig=plt.figure(figsize=(10,8))
fig.suptitle('star: '+arq[5:-4]+'   '+stpars)
minval=min(min(abon_i[0]), min(abon_ii[0]))-0.6
maxval=max(max(abon_i[0]), max(abon_ii[0]))+0.6

fig.add_axes([0.09, 0.55, 0.38, 0.37])
plt.axis([0,6,minval,maxval])
plt.plot(line1x_i, line1y_i, 'k-')
colorplot(chiex_i[0],abon_i[0], lbd_i[0])
#
plt.plot(rej_chiex,rej_abon, 'rx')
#
plt.text(0.04,0.03, '[FeI/H]='+'%6.4f'%abondi[0], fontsize=8, transform = plt.gca().transAxes)
plt.text(0.04,0.08, 'slope='+'%6.4f'%achiex_i[0], fontsize=8, transform = plt.gca().transAxes)
plt.figtext(0.28,0.50, r'$\chi_{exc}$ [eV]' ,fontdict={'fontsize':10}, horizontalalignment='center')
plt.figtext(0.025,0.73,"[X/H]",fontdict={'fontsize':10},rotation=90, verticalalignment='center')

fig.add_axes([0.58, 0.55, 0.38, 0.37])
plt.axis([0,6,minval,maxval])
plt.plot(line1x_ii, line1y_ii, 'k--')
colorplot(chiex_ii[0],abon_ii[0], lbd_ii[0])
plt.text(0.04,0.03, '[FeII/H]='+'%6.4f'%abondii[0], fontsize=8, transform = plt.gca().transAxes)
plt.figtext(0.77,0.50, r'$\chi_{exc}$ [eV]' ,fontdict={'fontsize':10}, horizontalalignment='center')
plt.figtext(0.518,0.73,"[X/H]",fontdict={'fontsize':10},rotation=90, verticalalignment='center')

fig.add_axes([0.09, 0.09, 0.38, 0.37])
plt.axis([-6.2 ,-4.6,minval,maxval])
plt.plot(line2x_i, line2y_i, 'k-')
colorplot(rdew_i[0],abon_i[0], lbd_i[0])
#
plt.plot(rej_ew,rej_abon, 'rx')
#
plt.text(0.04,0.03, '[FeI/H]='+'%6.4f'%abondi[0], fontsize=8, transform = plt.gca().transAxes)
plt.text(0.04,0.08, 'slope='+'%6.4f'%ardew_i[0], fontsize=8, transform = plt.gca().transAxes)
plt.figtext(0.28,0.04, 'red. EW [log($\AA$/$\lambda$)]' ,fontdict={'fontsize':10}, horizontalalignment='center')
plt.figtext(0.025,0.27,"[X/H]",fontdict={'fontsize':10},rotation=90, verticalalignment='center')


fig.add_axes([0.58, 0.09, 0.38, 0.37])
plt.axis([-6.2,-4.6,minval,maxval])
plt.plot(line2x_ii, line2y_ii, 'k--')
colorplot(rdrdew_ii[0],abon_ii[0], lbd_ii[0])
plt.text(0.04,0.03, '[FeII/H]='+'%6.4f'%abondii[0], fontsize=8, transform = plt.gca().transAxes)
plt.figtext(0.77,0.04, 'red. EW [log($\AA$/$\lambda$)]' ,fontdict={'fontsize':10}, horizontalalignment='center')
plt.figtext(0.518,0.27,"[X/H]",fontdict={'fontsize':10},rotation=90, verticalalignment='center')

plt.show()
print 'final met=%6.4f'%((abondii[0]+abondi[0])/2.)
