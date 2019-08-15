library(seqinr)

#goal, see if we can predict based on the distribution of ACGT in their respective genomes, if a piece of dengue fits dengue better than Zika!
#first lets make a transition probability matrix for both. These are the transitions of the Markov chain. Given the previous nucleotide, it can predict what future occurences may be.

#load data
zikaFasta <- read.fasta("zikaVirus.fasta")
dengueFasta <- read.fasta("dengue.fa")
zikaGenome <- zikaFasta[[1]]
dengueGenome <- dengueFasta[[1]]
zikaGenomeLength <- length(zikaGenome)
dengueGenomeLength <- length(dengueGenome)
#########################
# Step 1, fit genome of Zika to Markov chain model. Estimating its transition probability matrix
zikaDinuclOccurance <- table(zikaGenome[-zikaGenomeLength], zikaGenome[-1])
zikaDiNuclFreq <- zikaDinuclOccurance / rowSums(zikaDinuclOccurance)
print(zikaDinuclOccurance)
print(zikaDiNuclFreq)


# Transition Model for Dengue virus
dengueDinuclOccurance <- table(dengueGenome[-dengueGenomeLength], dengueGenome[-1])
dengueDinuclFreq <- dengueDinuclOccurance / rowSums(dengueDinuclOccurance)

print(dengueDinuclOccurance)
print(dengueDinuclFreq)

#########################

#next I"ll take 100 nucleotides from the dengue, and compare it to the probability in the MCM of the dengue virus to that of the Zika.
#Using log0-likelihood, we can determine that the sequence matches its own dna more closely than that of Zika's.

calculateBestMatchForSequence <- function(zikaTransMatrix, dengueTransMatrix, seq, order){ #calcs order
	likelihoodDecision <- 0
	seqLength <- length(seq)
	# Start at order + 1, as you want to compare the previous with the current nucleotide
	startNucleotide <- order + 1
	n = 0
	for (s in startNucleotide:seqLength) {
		# Take the range from
		n = n + 1
		#beginning, and ending
		beginKmer <- s - order
		endKmer <- s - 1
		# get position of nucleotides
		previousNucleotide <- seq[c(beginKmer:endKmer)]
		currentNucleotide <- seq[s]
		if (order > 1){
			previousNucleotide <- paste(previousNucleotide, collapse = ":")  #as previous nucleotide is dependant on order, length is longer in case of higher order
		} #get likeliood using formula log matrix prev/current. Div other matrix, same thing. Do for each
		likelihoodDecision <- likelihoodDecision + log(zikaTransMatrix[previousNucleotide, currentNucleotide] / dengueTransMatrix[previousNucleotide, currentNucleotide])
		}

	return(likelihoodDecision)
	}
#take 100 base pairs from position 101 to 200, could do other ofc
testSeq <- dengueGenome[101:200]
#now let's calculate it with order 1
decisionOrderOne <- calculateBestMatchForSequence(zikaDiNuclFreq, dengueDinuclFreq, testSeq, order =1)
print ('Order One')
print (decisionOrderOne)
#If log likelihood is negative, then it belongs to Dengue, else it belongs to Zika

if (decisionOrderOne >= 0){
	print("Sequence belongs to Zika")
} else if (decisionOrderOne == 0){
	print("Sequence belongs with a same degree to Zika and Dengue")
} else{
	print("Sequence belongs to Dengue")
}
####################3
# Now let's fit the Zica virus to a two second order Markov Chain and compare it to the simple MCM

getTrinucleotideFrequencies <- function(genome){ #second order
	genomeLength <- length(genome)
	oneBeforeTheLast <- genomeLength - 1
	twoBeforeLast <- genomeLength - 2

	firstnuclInTri <- factor(genome[c(1:twoBeforeLast)])
	secondnuclInTri <- factor(genome[c(2:oneBeforeTheLast)])
	thirdnuclInTri <- factor(genome[c(3:genomeLength)])
#	print('firstnuclInTri')
#	print(firstnuclInTri)
#	print('secondnuclInTri')
#	print(secondnuclInTri)
#	print('thirdnuclInTri')
#	print(thirdnuclInTri)
	trinuclOccurance <- table(firstnuclInTri:secondnuclInTri, thirdnuclInTri)
#	print('trinuclOccurance')
#	print(trinuclOccurance)
	trinuclFreq <- trinuclOccurance / rowSums(trinuclOccurance)
#	print('triNuclFreq')
#	print(trinuclFreq)
	return(trinuclFreq)
}

zikaTriNuclFreq <- getTrinucleotideFrequencies(zikaGenome)
dengueTriNuclFreq <- getTrinucleotideFrequencies(dengueGenome)
#now let's try second order
decisionOrderTwo <- calculateBestMatchForSequence(zikaTransMatrix = zikaTriNuclFreq, dengueTransMatrix =
dengueTriNuclFreq, seq = testSeq, order = 2)

print('zikaTriNuclFreq')
print(zikaTriNuclFreq)
print('dengueTriNuclFreq')
print(dengueTriNuclFreq)

print (decisionOrderTwo)


#same as before, is negative it's Dengue, else it's Zika
if (decisionOrderTwo >= 0){
	print("Sequence belongs to Zika")
}else if (decisionOrderTwo == 0){
	print("Sequence belongs with a same degree to Zika and Dengue")
}else{
	print("Sequence belongs to Dengue")
}

# Compare Values
sprintf("Difference between one and two order Malkov Chains: %f", {decisionOrderOne - decisionOrderTwo})


#for fun let's do third order, aka given three, get fourth
getQuatnucleotideFrequencies <- function(genome){
	genomeLength <- length(genome)
	oneBeforeTheLast <- genomeLength - 1
	twoBeforeLast <- genomeLength - 2
	thirdBeforeLast <- genomeLength - 3
	firstnuclInTri <- factor(genome[c(1:thirdBeforeLast)])
	secondnuclInTri <- factor(genome[c(2:twoBeforeLast)])
	thirdnuclInTri <- factor(genome[c(3:oneBeforeTheLast)])
	fourthNuclInTri <- factor(genome[c(4:genomeLength)])
	trinuclOccurance <- table(firstnuclInTri:secondnuclInTri:thirdnuclInTri, fourthNuclInTri)
	trinuclFreq <- trinuclOccurance / rowSums(trinuclOccurance)


#	print('trinuclOccurance')
#	print(trinuclOccurance)
#	print('trinuclFreq')
#	print(trinuclFreq)
	return(trinuclFreq)
}

zikaQuatNuclFreq <- getQuatnucleotideFrequencies(zikaGenome)
dengueQuatNuclFreq <- getQuatnucleotideFrequencies(dengueGenome)
decisionOrderThree <- calculateBestMatchForSequence(zikaTransMatrix = zikaQuatNuclFreq,
dengueTransMatrix = dengueQuatNuclFreq, seq = testSeq, order = 3)
print(decisionOrderThree)

#same as before, is negative it's Dengue, else it's Zika, this one should continue the trend
if (decisionOrderThree >= 0){
	print("Sequence belongs to Zika")
} else if (decisionOrderThree == 0){
	print("Sequence belongs with a same degree to Zika and Dengue")
} else{
	print("Sequence belongs to Dengue")
}

#yep belongs to Dengue

#now let's calculate the BIC, or bayesian information criterion, which should tell us which model is better for predicting on the Zika virus

##Calculating Best Order
#let's reload Zika
zikaFasta <- read.fasta("zikaVirus.fasta")
zikaGenome <- zikaFasta[[1]]


##For Order 0
n=length(zikaGenome)
par=3
c=count(zikaGenome,1)
#print(c)
p=count(zikaGenome,1,freq=T)
BIC=-2*sum(c*log(p)) + par*log(n)
print(BIC)

##For Order 1
n=length(zikaGenome)-1
par=12
a=count(zikaGenome,2)
print(a)
c = matrix(a, 4, 4, byrow=TRUE, dimnames = list(c("A","C","G","T"),c("A","C","G","T")))
#print(c)
p=c[,]/(c[,1]+c[,2]+c[,3]+c[,4])
BIC=-2*sum(c*log(p)) + par*log(n)
print(BIC)


##For Order 2
n=length(zikaGenome)-2
print (n)
par=48
a=count(zikaGenome,3)
c = matrix(a, 16, 4, byrow=TRUE, dimnames = list(c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"),c("A","C","G","T")))
#print(c)
p=c[,]/(c[,1]+c[,2]+c[,3]+c[,4])
BIC=-2*sum(c*log(p)) + par*log(n)
print(BIC)

#for fun can do order 3
##for Order k (3)
n=length(zikaGenome)-3
par=196#
a=count(zikaGenome,4)
#had python print out the list for the next one rather than type it by hand.
c = matrix(a, 64, 4, byrow=TRUE, dimnames = list(c('AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CAT','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAA','TAC','TAG','TAT','TCA','TCC','TCG','TCT','TGA','TGC','TGG','TGT','TTA','TTC','TTG','TTT'),c("A","C","G","T")))
#print(c)
p=c[,]/(c[,1]+c[,2]+c[,3]+c[,4])
BIC=-2*sum(c*log(p)) + par*log(n)
print(BIC)

##for Order k (4)
#n=length(zikaGenome)-4
#par=584
#c=table(LargeTable)
#print(c)
#p=c[,]/(c[,1]+c[,2]+c[,3]+c[,4])
#BIC=-2*sum(c*log(p)) + par*log(n)
#print(BIC)


#using first order, the decision criterion for dengue was 1.496 over that of Zika.
#2nd and 3rd order are also calculated, and gave 2.772 and 5.422 respectively.


#Bayesian information criterion (BIC) can be used to select the higher models. A bic of 29115 was found for order 1, which is better than that of order 0 of 29751. Also better than that of Order 2 which is 29244.2

#for fun showed order 3, which seems to be 30384, I started on order four but there's no point.

