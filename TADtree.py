import numpy as np, scipy.stats as ss, os, sys
#----------------------------------------------------------------------------------------#
#							  		READ CONTROL FILE	
#----------------------------------------------------------------------------------------#
chrs = [] 
paths = []
N = []
S = None
p = None
q = None
M = None
gamma = None
output_directory = None

for l in open(sys.argv[1]).read().split('\n'):
	l = l.split()
	if len(l) > 2:
		if l[0] == 'S' and l[1] == '=': S = int(l[2])
		if l[0] == 'p' and l[1] == '=': p = int(l[2])
		if l[0] == 'q' and l[1] == '=': q = int(l[2])
		if l[0] == 'M' and l[1] == '=': M = int(l[2])
		if l[0] == 'gamma' and l[1] == '=': gamma = int(l[2])
		if l[0] == 'output_directory' and l[1] == '=': output_directory = l[2]
		if l[0] == 'contact_map_path' and l[1] == '=': paths = [x for x in l[2].split(',')]
		if l[0] == 'contact_map_name' and l[1] == '=': chrs = [x for x in l[2].split(',')]
		if l[0] == 'N' and l[1] == '=': N = [int(x) for x in l[2].split(',')]

if S == None: print 'Havent specified parameter: "S"'
if p == None: print 'Havent specified parameter: "p"'
if q == None: print 'Havent specified parameter: "q"'
if M == None: print 'Havent specified parameter: "M"'
if gamma == None: print 'Havent specified parameter: "gamma"'
if output_directory == None: print 'Havent specified parameter: "output_directory"'

if len(chrs) != len(paths) or len(N) != len(paths): print 'Number of contact map paths, contact map names, and N values is inconsistent'
if len(paths) == 0: print 'No contact map paths specified'

#----------------------------------------------------------------------------------------#
#							   LOAD CONTACTS AND BACKGROUND
#----------------------------------------------------------------------------------------#
# load data
print 'Loading data'
mats = {chrs[i] : np.loadtxt(paths[i]) for i in range(len(paths))}
height = S
backbins = []
for chr in chrs:
	for i in range(mats[chr].shape[0]-height):
		backbins += [mats[chr][i,i:i+height]]
backgrnd = np.mean(backbins,axis=0)

def normalize(mat, background):
	out = np.array(mat)
	n = len(background)
	for i in range(mat.shape[0]):
		for j in range(max([0,i-n+1]),min([i+n,mat.shape[0]])):
			if mat[i,j] > 0:
				out[i,j] = np.log(mat[i,j] / background[int(np.abs(j-i))])
			else: out[i,j] = 0
	return out
#----------------------------------------------------------------------------------------#
#									PRECOMPUTE SCORES
#----------------------------------------------------------------------------------------#
def betadelta(chr,i,j):
	n = j-i
	x,y = [],[]
	for k in range(n):
		for l in range(k+1,n):
			x += [float(l-k)]
			y += [mats[chr][i+k,i+l] / backgrnd[l-k]]
	x = np.array(x)
	y = np.array(y)
	delta, beta, r_value, p_value, std_err = ss.linregress(x,y)

	if delta < 0:  fit = np.inf
	else:
		fit = 0
		for k in range(n-1):
			for l in range(k+1,n):
				fit += ((float(l-k)*delta+beta)*backgrnd[l-k] - mats[chr][i+k,i+l])**2
	return beta,delta,fit

#----------------------------------------------------------------------------------------#
tadscores = {}
bakscores = {}
chrdeltas = {}
chrbetas = {}

print 'Precomputing paramters for ...'
for chr in chrs:
	print chr
	n = mats[chr].shape[0]
	smat = np.zeros((n,n))
	gmat = np.zeros((n,n))
	bmat = np.zeros((n,n))
	for i in range(n-2):
		for j in range(i+3,i+np.min([n-i, height])):
			beta,delta,fit = betadelta(chr,i,j)
			smat[i,j] = fit
			smat[j,i] = fit
			gmat[i,j] = delta
			bmat[i,j] = beta
	tadscores.update({chr:smat})
	chrdeltas.update({chr:gmat})
	chrbetas.update({chr:bmat})

print '\nPrecomputing background scores for ...'
for chr in chrs:
	print chr
	n = mats[chr].shape[0]
	smat = np.zeros((n,n))
	for i in range(n-2):
		for j in range(i+2,i+np.min([n-i, height])):
			fit = 0
			for k in range(j-i-1):
				for l in range(k+1,j-i):
					fit += (mats[chr][i+k,i+l]-backgrnd[l-k])**2
			smat[i,j] = fit 
			smat[j,i] = fit
	bakscores.update({chr:smat})
#----------------------------------------------------------------------------------------#
#				BOUNDED HIERARCHICAL WEIGHTED INTERVAL SCHEDULING
#----------------------------------------------------------------------------------------#	
def buildtrees(mat,smat,gmat,bmat,bakmat,t_lim,height,min_size):
	L = smat.shape[0]
	if height > L: height = L
	score = np.zeros((L,L,t_lim))
	traceback = np.zeros((L,L,t_lim),dtype=int)
	local_partitions = {}
	
	for n in range(min_size,height):
		for i in range(L - n):			
			j = i + n
			if smat[i,j] > 0: continue
			if gmat[i,j] < 0:
				score[i,j,:] = np.inf
			else: 
				score[i,j,0] = smat[i,j]
							
				local_t_lim = np.min([t_lim,n-min_size+1])
				local_score = np.zeros((n+1,local_t_lim))
				local_traceback_k = np.zeros((n+1,local_t_lim),dtype=int)
				local_traceback_t = np.zeros((n+1,local_t_lim),dtype=int)
			
				for k in range(min_size,n+1):
					for t in range(1,np.min([t_lim,k-min_size+1])):
						options = np.zeros((k+1,t))
						options[0,0] = local_score[k-1,t]
						for l in range(min_size,k+1):
							for tt in range(np.min([l-min_size,t])):				
								olddelta = gmat[i,j]
								oldbeta = bmat[i,j]
								if olddelta > gmat[i+k-l,i+k]:
									options[l,tt] = np.inf
								else:
									oldscore = 0
									for z in range(l-1):
										for w in range(z+1,l):
											oldscore += ((olddelta*float(w-z)+oldbeta)*backgrnd[w-z] - mat[i+k-l+z,i+k-l+w])**2
							
									options[l,tt] = local_score[k-l,t-tt-1] + score[i+k-l,i+k,tt] - (oldscore - bakmat[i+k-l,i+k])
						best = np.argmin(options)
						local_score[k,t] = np.min(options)
						if best == 0:
							local_traceback_k[k,t] = -1
							local_traceback_t[k,t] = t
						else:
							best_l = best / t
							best_t = best % t
							local_traceback_k[k,t] = best_l
							local_traceback_t[k,t] = best_t
				
				for t in range(1,local_t_lim):
					score[i,j,t] = local_score[n,t] + smat[i,j]
					intervals = []
					pos = n
					tt = t
					while True:
						if pos <= min_size or local_traceback_k[pos,tt] == 0: break
						if local_traceback_k[pos,tt] == -1:
							pos = pos-1
						else:
							newpos = pos-local_traceback_k[pos,tt]
							intervals += [[i+newpos,i+pos,local_traceback_t[pos,tt]]]
							newtt = tt-local_traceback_t[pos,tt]-1
							tt = newtt
							pos = newpos
					local_partitions.update({(i,j,t):np.array(intervals)})
	local_parts_array = np.zeros((L,height,t_lim,t_lim,3),dtype=int)
	for triple in local_partitions:
		i,j,t = triple
		l = len(local_partitions[triple])
		if l > 0:
			local_parts_array[i,j-i,t,:l,:] = local_partitions[triple]
		
	return local_parts_array,score	

def getforest(score,height,T_lim,t_lim,min_size):
	L = mat.shape[1]
	totalscore = np.zeros((L,T_lim))
	traceback_k = np.zeros((L,T_lim),dtype=int)
	traceback_t = np.zeros((L,T_lim),dtype=int)
	for i in range(min_size,L):
		for t in range(1,T_lim):
			options = np.zeros((i+1,t))
			options[0,0] = totalscore[i-1,t]
			for k in range(min_size,np.min([height,i])):		
				for tt in range(np.min([t,t_lim])):		
					options[k,tt] = score[i-k,i,tt] + totalscore[i-k,t-tt-1]

			totalscore[i,t] = np.min(options)
			best = np.argmin(options)
			if best == 0:
				traceback_k[i,t] = -1
				traceback_t[i,t] = t
			else:
				best_k = best / t
				best_t = best % t			
				traceback_k[i,t] = best_k
				traceback_t[i,t] = best_t
	return totalscore,traceback_k, traceback_t

def foresttb(totalscore,traceback_k, traceback_t,start_t):
	L = totalscore.shape[0]
	trees = []
	pos = L-1
	tt = start_t
	while True:
		if pos <= min_size or tt == 0: break
		if traceback_k[pos,tt] == -1:
			pos = pos - 1
		else:
			newpos = pos - traceback_k[pos,tt]
			trees += [(newpos,pos,traceback_t[pos,tt])]
			tt = tt - traceback_t[pos,tt] - 1
			pos = newpos
	return trees

def all_intervals(local_parts_array,i,j,t):
	intervals = [(i,j,t)]
	if t > 0:
		for p in local_parts_array[i,j-i,t]:
			if np.sum(p) > 0:
				intervals += all_intervals(local_parts_array,p[0],p[1],p[2])
	return intervals				
#----------------------------------------------------------------------------------------#
print '\nRunning dynamic program for ...'
for i,chr in enumerate(chrs):
	print chr
	min_size = 2
	t_lim = M
	T_lim = N[i]

	mat = mats[chr]
	gmat = chrdeltas[chr]
	bmat = chrbetas[chr]
	bakmat = bakscores[chr]
	smat = tadscores[chr] - bakscores[chr]
	
	###############
	short = p; long = 1; steps=q
	bi = []
	for i in range(long,mat.shape[0]-long):
		b = 0
		for s in range(1,steps):
			a1 = np.sum(mat[i-long*s:i-long*(s-1),i-short:i])
			b1 = np.sum(mat[i-long*s:i-long*(s-1),i:i+short])
			a2 = np.sum(mat[i+long*(s-1):i+long*s,i-short:i])
			b2 = np.sum(mat[i+long*(s-1):i+long*s,i:i+short])	
			b += np.abs(a1-b1) + np.abs(a2-b2)	
		bi += [b]
	bi = [0]*long + bi + [0]*long
	bi = (np.array(bi) - np.mean(bi)) / np.std(bi) 

	for i in range(mat.shape[0]):
		for j in range(mat.shape[1]):
			if bi[i] < 0 or bi[j] < 0: smat[i,j] = np.inf
			else: smat[i,j] = smat[i,j] - gamma*(bi[i]+bi[j])
			
	
	###############
	print 'Building TADtrees for ' + chr
	local_parts_array,score = buildtrees(mat,smat,gmat,bmat,bakmat,t_lim,height,min_size)
	print 'Assembling TADforest for ' + chr
	totalscore,traceback_k,traceback_t = getforest(score,height,T_lim,t_lim,min_size)

	if not os.path.exists(output_directory + '/' + chr):
		os.system('mkdir ' + output_directory + '/' + chr)
	duplicates_out = ['\t'.join(['name','proportion_duplicates'])]

	for start_t in range(1,T_lim):
		trees = foresttb(totalscore,traceback_k, traceback_t,start_t)
		allints = []

		for t in trees: 
			allints += all_intervals(local_parts_array,t[0],t[1],t[2])
	
		allints = sorted(allints, key = lambda x: np.abs(x[1]-x[0]))
	
		final_ints = []
		for t in allints:
			duplicate = False
			for tt in final_ints:
				if (np.abs(t[0] - tt[0]) < 2 and np.abs(t[1]-tt[1]) < 2):
					duplicate = True
			if not duplicate:
				final_ints += [[t[0],t[1]]]
				
					
		out = ['\t'.join(['chr','start','end'])]
		for t in sorted(final_ints, key=lambda x: x[0]):
			out += ['\t'.join([chr,repr(t[0]),repr(t[1])])]		
		fname = output_directory + '/' + chr + '/N' + repr(start_t) + '.txt'
		open(fname,'w').write('\n'.join(out))
		
		duplicates_out += ['\t'.join(['N'+repr(start_t), repr(1 - float(len(final_ints)) / len(allints))])]
		
	open(output_directory + '/' + chr + '/proportion_duplicates.txt','w').write('\n'.join(duplicates_out))