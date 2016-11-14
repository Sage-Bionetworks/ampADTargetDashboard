#push target list to synapse, share with Kenny

#foo = read.csv('targetsV1.csv',stringsAsFactors=F,header = F)

#library(utilityFunctions)

#bar = utilityFunctions::convertHgncToEnsembl(foo)

#map to ensembl id

require(synapseClient)
synapseLogin()
foo = synGet('syn7525109')
bar = read.csv(foo@filePath,stringsAsFactors=F)
View(bar)

#figure out which network to pull down

netObj = synGet('syn5909604')
load(netObj@filePath)
sum(bicNetworks$rankConsensus$network)

#extract local neighborhoods
library(dplyr)
net1 = bicNetworks$rankConsensus$network %>% as.matrix
net1 = net1 + t(net1)


spikeIn = colnames(net1)%in% bar$ENSG
neighb1 = net1%*%spikeIn !=0

#create custom edge list 
keepGene = colnames(net1)[which(neighb1==1)]
net1 = net1[keepGene,keepGene]

edgeList = metanetwork::rankedEdgeList(net1,symmetric=TRUE)
edgeList2 = edgeList
edgeList2$var1 = colnames(net1)[edgeList$var1]
edgeList2$var2 = colnames(net1)[edgeList$var2]
edgeList2 = dplyr::select(edgeList2,var1,var2)

write.csv(edgeList2,file='tempNetWoT.csv',quote=F,row.names=F)

tempNetWoTObj = File('tempNetWoT.csv',parentId='syn7525109')


#share with Kenny


#pull down expression for local neighborhoods

#pull down clinical assessment for patients

#share with kenny appropriately