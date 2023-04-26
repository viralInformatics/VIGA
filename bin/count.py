import sys
indir=sys.argv[1]
consensus=sys.argv[2]

pileup = open(indir, "r")
out=open(consensus+'.txt', "w")

out.write('Chr\tPosition\tA_num\tT_num\tC_num\tG_num\tTrue\tRef\n')
for line in pileup:
	linesplit = line.split("\t")
	chrom = linesplit[0]
	position = linesplit[1]
	ref = linesplit[2]
	A =  str(linesplit[4].count("A") + linesplit[4].count("a"))
	T =  str(linesplit[4].count("T") + linesplit[4].count("t"))
	C =  str(linesplit[4].count("C") + linesplit[4].count("c"))
	G =  str(linesplit[4].count("G") + linesplit[4].count("g"))

	if ref == "A":
		A =  str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "T":
		T =  str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "C":
		C =  str(linesplit[4].count(".") + linesplit[4].count(","))
	elif ref == "G":
		G =  str(linesplit[4].count(".") + linesplit[4].count(","))
	else:
		continue
	
	if max(int(A),int(T),int(C),int(G))==0:
		output = "\t".join([chrom, position, A, T, C, G, 'N',ref])
	else:
		if max(int(A),int(T),int(C),int(G))==int(A):
			output = "\t".join([chrom, position, A, T, C, G, 'A',ref])
		elif max(int(A),int(T),int(C),int(G))==int(T):
			output = "\t".join([chrom, position, A, T, C, G, 'T',ref])
		elif max(int(A),int(T),int(C),int(G))==int(C):
			output = "\t".join([chrom, position, A, T, C, G, 'C',ref])
		elif max(int(A),int(T),int(C),int(G))==int(G):
			output = "\t".join([chrom, position, A, T, C, G, 'G',ref])
	out.write(output+'\n')

fa=open(consensus+'.txt','r')
faout=open(consensus,'w')
faout.write('>consensus_fasta\n')
next(fa)
for line in fa:
	linesplit = line.split("\t")
	con=linesplit[6]
	faout.write(con)
faout.write('\n')
pileup.close()
out.close()
fa.close()
faout.close()
