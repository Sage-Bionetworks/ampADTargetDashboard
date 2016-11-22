#push noam's networks to Synapse

require(synapseClient)
synapseLogin()

system('unzip reurgentdatastillneededformssmcleanup.zip')
foo = system(command = 'ls',intern = TRUE)

keepFile = grep('gene',foo)

bar = lapply(foo[keepFile],function(x) return(read.delim(x,stringsAsFactors=F,header=F,sep='\t')))
bar[[1]]
foobar = do.call(rbind,bar)
View(foobar)

colnames(foobar)[1:2] = c('var1','var2')
write.csv(foobar[,1:2],file='mssmNetwork.csv',quote=F,row.names=F)

bax = File('mssmNetwork.csv',parentId='syn7525089')

permLink =githubr::getPermlink('Sage-Bionetworks/WallOfTargets','uploadMSSMNetworks.R')

bax = synStore(bax,used = permLink)
