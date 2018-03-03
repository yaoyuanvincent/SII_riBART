###This code contains the functions needed to run BART.

#Author: Yaoyuan Vincent, Tan
#Date created: 31 May 2016
#Date updated: 20 Jul 2016

growstep=function(variable,covmat,value,xindex,varindex,depth,outcome){

	#set all outputs to the default values first
	variable_grow=variable
	value_grow=value
	xindex_grow=xindex
	varindex_grow=varindex
	depth_grow=depth
	num_term=1
	num_pred=ncol(covmat)
	num_value=0
	num_obs_in_pick_node=0
	num_obs_in_left_child=0
	num_obs_in_right_child=0
	out_in_left_child=NULL
	out_in_right_child=NULL
	out_in_pick_node=NULL
	depth_pick=1
	num_2g_int_nodes=1

	#check whether it is a single node
	if (length(value)==1){
		col_available=NULL
		for (i in 1:ncol(covmat)){
			if(length(unique(covmat[,i]))>1){
				col_available=c(col_available,i)
			}
		}
		if(length(col_available)>0){
			if(length(col_available)==1){
				pick_col=col_available
			} else {
				pick_col=sample(col_available,1)
			}
			raw_var=covmat[,pick_col]
			if(length(unique(raw_var))>2){
				sample_value=sort(unique(raw_var))
				sample_value=sample_value[-1]
				sample_value=sample_value[-length(sample_value)]
				if(length(sample_value)==1){
					pick_value=sample_value
				} else {
					pick_value=sample(sample_value,1)
				}			
				value_grow=c(pick_value,NA,NA)
				xindex_grow=c(pick_col,NA,NA)
				depth_grow=c(1,2,2)
				variable_grow=cbind(ifelse(raw_var<pick_value,1,0),1-ifelse(raw_var<pick_value,1,0))
				varindex_grow=c(2,3)
				num_pred=length(col_available)
				num_value=length(unique(raw_var))
				num_obs_in_pick_node=length(raw_var)
				num_obs_in_left_child=sum(variable_grow[,1])
				num_obs_in_right_child=sum(variable_grow[,2])
				out_in_left_child=outcome[variable_grow[,1]==1]
				out_in_right_child=outcome[variable_grow[,2]==1]
				out_in_pick_node=outcome
			}
		}
	} else {
		pick_node=sample(ncol(variable),1)
		new_covmat=covmat[variable[,pick_node]==1,]
		if(!is.vector(new_covmat)){
			if(nrow(new_covmat)>2){
				col_available=NULL
				for (i in 1:ncol(new_covmat)){
					if(length(unique(new_covmat[,i]))>1){
						col_available=c(col_available,i)
					}
				}
				if(length(col_available)>0){
					if(length(col_available)==1){
						pick_col=col_available
					} else {
						pick_col=sample(col_available,1)
					}
					raw_var=new_covmat[,pick_col]
					if(length(unique(raw_var))>2){
						sample_value=sort(unique(raw_var))
						sample_value=sample_value[-1]
						sample_value=sample_value[-length(sample_value)]
						if(length(sample_value)==1){
						pick_value=sample_value
						} else {
						pick_value=sample(sample_value,1)
						}
						if (varindex[pick_node]!=length(value)){
							value_grow=c(value[1:(varindex[pick_node]-1)],pick_value,NA,NA,value[(varindex[pick_node]+1):length(value)])
							xindex_grow=c(xindex[1:(varindex[pick_node]-1)],pick_col,NA,NA,xindex[(varindex[pick_node]+1):length(xindex)])
							depth_grow=c(depth[1:(varindex[pick_node])],depth[varindex[pick_node]]+1,depth[varindex[pick_node]]+1,depth[(varindex[pick_node]+1):length(depth)])
						} else {
							value_grow=c(value[1:(varindex[pick_node]-1)],pick_value,NA,NA)
							xindex_grow=c(xindex[1:(varindex[pick_node]-1)],pick_col,NA,NA)
							depth_grow=c(depth,depth[length(value)]+1,depth[length(value)]+1)
						}
						if (pick_node==1){
							variable_grow=cbind(variable[,1]*ifelse(covmat[,pick_col]<pick_value,1,0),variable[,1]*(1-ifelse(covmat[,pick_col]<pick_value,1,0)),variable[,2:ncol(variable)])
							varindex_grow=c(varindex[1]+1,varindex[1]+2,varindex[2:length(varindex)]+2)
						} else if(pick_node==ncol(variable)){
							variable_grow=cbind(variable[,1:(ncol(variable)-1)],variable[,pick_node]*ifelse(covmat[,pick_col]<pick_value,1,0),variable[,pick_node]*(1-ifelse(covmat[,pick_col]<pick_value,1,0)))
							varindex_grow=c(varindex[1:(pick_node-1)],varindex[pick_node]+1,varindex[pick_node]+2)
						} else {
							variable_grow=cbind(variable[,1:(pick_node-1)],variable[,pick_node]*ifelse(covmat[,pick_col]<pick_value,1,0),variable[,pick_node]*(1-ifelse(covmat[,pick_col]<pick_value,1,0)),variable[,(pick_node+1):ncol(variable)])
							varindex_grow=c(varindex[1:(pick_node-1)],varindex[pick_node]+1,varindex[pick_node]+2,varindex[(pick_node+1):length(varindex)]+2)
						}
						num_2g_int_nodes=0
						for(i in 1:(length(value_grow)-2)){
							if(!is.na(value_grow[i])&is.na(value_grow[i+1])&is.na(value_grow[i+2])){
								num_2g_int_nodes=num_2g_int_nodes+1
							}
						}
						left_child=variable[,pick_node]*ifelse(covmat[,pick_col]<pick_value,1,0)
						right_child=variable[,pick_node]*(1-ifelse(covmat[,pick_col]<pick_value,1,0))
						num_term=ncol(variable)
						out_in_pick_node=outcome[variable[,pick_node]==1]
						num_pred=length(col_available)
						num_value=length(unique(raw_var))
						num_obs_in_pick_node=nrow(new_covmat)
						num_obs_in_left_child=sum(left_child)
						num_obs_in_right_child=sum(right_child)
						out_in_left_child=outcome[left_child==1]
						out_in_right_child=outcome[right_child==1]
						depth_pick=depth[varindex[pick_node]]		
					}
				}
			}
		}
	}
	return(list(variable=variable_grow,varindex=varindex_grow,value=value_grow,xindex=xindex_grow,num_term=num_term
		,num_pred=num_pred,num_value=num_value,num_obs_in_pick_node=num_obs_in_pick_node
		,num_obs_in_left_child=num_obs_in_left_child,num_obs_in_right_child=num_obs_in_right_child,out_in_left_child=out_in_left_child
		,out_in_right_child=out_in_right_child,out_in_pick_node=out_in_pick_node
		,depth_pick=depth_pick,num_2g_int_nodes=num_2g_int_nodes,depth=depth_grow))
}

MHgrow=function(prob_prune,prob_grow,alpha,beta,sigma,sigma_mu,num_term,num_pred,num_value,num_2g_int_nodes,num_obs_in_pick_node,num_obs_in_left_child,num_obs_in_right_child,out_in_left_child,out_in_right_child,out_in_pick_node,depth_pick){
	if(!is.null(out_in_left_child)){
		#draw from uniform distirbution
		logu=log(runif(1))
		#According to Dr. Kapelner's paper, r can be decomposed into 3 ratios
		#Transition ratio
		ltr=log(prob_prune)-log(prob_grow)+log(num_term)+log(num_pred)+log(num_value)-log(num_2g_int_nodes)
		#Likelihood ratio
		llr=0.5*2*log(sigma)+log(sigma^2+num_obs_in_pick_node*sigma_mu^2)-log(sigma^2+num_obs_in_left_child*sigma_mu^2)-log(sigma^2+num_obs_in_right_child*sigma_mu^2)+(sigma_mu^2/(2*sigma^2))*(((sum(out_in_left_child))^2)/(sigma^2+num_obs_in_left_child*sigma_mu^2)+((sum(out_in_right_child))^2)/(sigma^2+num_obs_in_right_child*sigma_mu^2)-((sum(out_in_pick_node))^2)/(sigma^2+num_obs_in_pick_node*sigma_mu^2))
		#Tree_structure_ratio
		ltsr=log(alpha)+2*log(1-alpha/((2+depth_pick)^beta))-log((1+depth_pick)^beta-alpha)-log(num_pred)-log(num_value)		
		logr=ltr+llr+ltsr
		tmpaccept=ifelse(logu<logr,1,0)
		accept=ifelse(is.na(tmpaccept),0,tmpaccept)
	} else {
		accept=0
		logu=NA
		ltr=NA
		llr=NA
		ltsr=NA
		logr=NA
	}
	return(list(accept=accept,lu=logu,ltr=ltr,llr=llr,ltsr=ltsr,lr=logr))
}

prunestep = function(variable,covmat,value,xindex,varindex,depth,outcome){
	#set all outputs to the default values first
	variable_prune=variable
	value_prune=value
	xindex_prune=xindex
	varindex_prune=varindex
	depth_prune=depth
	num_term=1
	num_pred=ncol(covmat)
	num_value=0
	num_obs_in_pick_node=0
	num_obs_in_left_child=0
	num_obs_in_right_child=0
	out_in_left_child=NULL
	out_in_right_child=NULL
	out_in_pick_node=NULL
	depth_pick=1
	num_2g_int_nodes=1
	#Do only if it is not a single node
	if (length(value)>1){
		if (length(value)==3){
			# because we know what this is just 1 internal node and 2 leaf nodes.
			variable_prune=rep(1,length(outcome))
			value_prune=NA
			xindex_prune=NA
			varindex_prune=1
			depth_prune=1
			#number of terminal nodes in new tree
			num_term=1
			#number of columns left in the new pruned tree to split on
			num_pred=ncol(covmat)
			#number of unique elements in the chosen variable that we can use to split
			num_value=length(unique(covmat[,xindex[1]]))
			num_obs_in_pick_node=length(outcome)
			num_obs_in_left_child=sum(variable[,1])
			num_obs_in_right_child=sum(variable[,2])
			out_in_left_child=outcome[variable[,1]==1]
			out_in_right_child=outcome[variable[,2]==1]
			out_in_pick_node=outcome
			depth_pick=1
			#number of singly nodes in the old tree
			num_2g_int_nodes=1
		} else {
			#Find the internal node to prune.
			avail_index=NULL
			num_2g_int_nodes=0
			for(i in 1:(length(value)-2)){
				if(!is.na(value[i])&is.na(value[i+1])&is.na(value[i+2])){
					num_2g_int_nodes=num_2g_int_nodes+1
					avail_index=c(avail_index,i)
				}
			}
			if(length(avail_index)==1){
				pick_index=avail_index
			} else {
				pick_index=sample(avail_index,1)
			}
			idxvar=1:ncol(variable)
			pidx=idxvar[varindex==(pick_index+1)]
			if (idxvar[varindex==(pick_index+1)]==1) {
				variable_prune=cbind(variable[,1]+variable[,2],variable[,3:ncol(variable)])
				varindex_prune=c(varindex[1]-1,varindex[3:length(varindex)]-2)
			} else if (idxvar[varindex==(pick_index+1)]==(length(idxvar)-1)) {
				variable_prune=cbind(variable[,1:(ncol(variable)-2)],variable[,(ncol(variable)-1)]+variable[,ncol(variable)])
				varindex_prune=c(varindex[1:(length(varindex)-2)],varindex[length(varindex)-1]-1)
			} else {
				variable_prune=cbind(variable[,1:(pidx-1)],variable[,varindex==(pick_index+1)]+variable[,varindex==(pick_index+2)],variable[,(pidx+2):ncol(variable)])
				varindex_prune=c(varindex[1:(pidx-1)],varindex[pidx]-1,varindex[(pidx+2):length(varindex)]-2)
			}
			num_term=ncol(variable_prune)
			new_covmat=as.matrix(covmat[(variable[,varindex==(pick_index+1)]+variable[,varindex==(pick_index+2)])==1,])
			col_available=NULL
			for (i in 1:ncol(new_covmat)){
				if(length(unique(new_covmat[,i]))>1){
					col_available=c(col_available,i)
				}
			}
			num_pred=length(col_available)
			num_value=length(unique(new_covmat[,xindex[pick_index]]))
			num_obs_in_pick_node=sum(variable[,varindex==(pick_index+1)]+variable[,varindex==(pick_index+2)])
			num_obs_in_left_child=sum(variable[,varindex==(pick_index+1)])
			num_obs_in_right_child=sum(variable[,varindex==(pick_index+2)])
			out_in_left_child=outcome[variable[,varindex==(pick_index+1)]==1]
			out_in_right_child=outcome[variable[,varindex==(pick_index+2)]==1]
			out_in_pick_node=outcome[(variable[,varindex==(pick_index+1)]+variable[,varindex==(pick_index+2)])==1]
			depth_pick=depth[pick_index]
			if(length(value)==(pick_index+2)){
				value_prune=c(value[1:(pick_index-1)],NA)
				xindex_prune=c(xindex[1:(pick_index-1)],NA)
				depth_prune=depth[1:pick_index]
			} else {
				value_prune=c(value[1:(pick_index-1)],NA,value[(pick_index+3):length(value)])
				xindex_prune=c(xindex[1:(pick_index-1)],NA,xindex[(pick_index+3):length(xindex)])
				depth_prune=c(depth[1:pick_index],depth[(pick_index+3):length(depth)])
			}
		}
	}
	return(list(variable=variable_prune,varindex=varindex_prune,value=value_prune,xindex=xindex_prune,num_term=num_term
		,num_pred=num_pred,num_value=num_value,num_obs_in_pick_node=num_obs_in_pick_node
		,num_obs_in_left_child=num_obs_in_left_child,num_obs_in_right_child=num_obs_in_right_child,out_in_left_child=out_in_left_child
		,out_in_right_child=out_in_right_child,out_in_pick_node=out_in_pick_node
		,depth_pick=depth_pick,num_2g_int_nodes=num_2g_int_nodes,depth=depth_prune))
}

MHprune=function(prob_prune,prob_grow,alpha,beta,sigma,sigma_mu,num_term,num_pred,num_value,num_2g_int_nodes,num_obs_in_pick_node,num_obs_in_left_child,num_obs_in_right_child,out_in_left_child,out_in_right_child,out_in_pick_node,depth_pick){
	if(!is.null(out_in_left_child)){
		#draw from uniform distirbution
		logu=log(runif(1))
		#According to Dr. Kapelner's paper, r can be decomposed into 3 ratios
		#Transition ratio
		ltr=log(prob_grow)-log(prob_prune)+log(num_2g_int_nodes)-log(num_term)-log(num_pred)-log(num_value)	
		#Likelihood ratio
		llr=0.5*2*log(sigma)+log(sigma^2+num_obs_in_pick_node*sigma_mu^2)-log(sigma^2+num_obs_in_left_child*sigma_mu^2)-log(sigma^2+num_obs_in_right_child*sigma_mu^2)+(sigma_mu^2/(2*sigma^2))*(((sum(out_in_left_child))^2)/(sigma^2+num_obs_in_left_child*sigma_mu^2)+((sum(out_in_right_child))^2)/(sigma^2+num_obs_in_right_child*sigma_mu^2)-((sum(out_in_pick_node))^2)/(sigma^2+num_obs_in_pick_node*sigma_mu^2))
		#Tree_structure_ratio
		ltsr=log(alpha)+2*log(1-alpha/((2+depth_pick)^beta))-log((1+depth_pick)^beta-alpha)-log(num_pred)-log(num_value)
		logr=ltr-llr-ltsr
		tmpaccept=ifelse(logu<logr,1,0)
		accept=ifelse(is.na(tmpaccept),0,tmpaccept)
	} else {
		accept=0
		logu=NA
		ltr=NA
		llr=NA
		ltsr=NA
		logr=NA
	}
	return(list(accept=accept,lu=logu,ltr=ltr,llr=llr,ltsr=ltsr,lr=logr))
}

sig_draw=function(y,predy,nu,lambda){
	#needed to draw inverse gamma
	library(pscl)
	return(sqrt(rigamma(1,(length(y)+nu)/2,(nu*lambda+sum((y-predy)^2))/2)))

}

bart1draw=function(transy,mu_ij,initsigma,sigma_mu,X,value,xindex,depth,varindex,variable,m,alpha,beta,probgrow,probprune,beta0,beta1){
	
	tmpvalue=list()
	tmpxindex=list()
	tmpdepth=list()
	tmpvarindex=list()
	tmpvariable=list()
	tmpmu=matrix(0,nrow=length(transy),ncol=m)
	mu_ij_list=list()
	
	ioutcome=transy-rowSums(mu_ij[,-1])
	ivalue=value[[1]]
	ixindex=xindex[[1]]
	idepth=depth[[1]]
	ivarindex=varindex[[1]]
	ivariable=variable[[1]]
	
	proposal_type=sample(2,1)

	if(proposal_type==1){
		tmpgrowstep=growstep(ivariable,X,ivalue,ixindex,ivarindex,idepth,ioutcome)
		tmpgrowMH=MHgrow(probgrow,probprune,alpha,beta,initsigma,sigma_mu,tmpgrowstep$num_term,tmpgrowstep$num_pred,tmpgrowstep$num_value,tmpgrowstep$num_2g_int_nodes,tmpgrowstep$num_obs_in_pick_node,tmpgrowstep$num_obs_in_left_child,tmpgrowstep$num_obs_in_right_child,tmpgrowstep$out_in_left_child,tmpgrowstep$out_in_right_child,tmpgrowstep$out_in_pick_node,tmpgrowstep$depth_pick)
	
		if(tmpgrowMH$accept==1){
			tmpvalue[[1]]=tmpgrowstep$value
			tmpxindex[[1]]=tmpgrowstep$xindex
			tmpdepth[[1]]=tmpgrowstep$depth
			tmpvarindex[[1]]=tmpgrowstep$varindex
			tmpvariable[[1]]=tmpgrowstep$variable	
		}else{
			tmpvalue[[1]]=ivalue
			tmpxindex[[1]]=ixindex
			tmpdepth[[1]]=idepth
			tmpvarindex[[1]]=ivarindex
			tmpvariable[[1]]=ivariable
		}
	} else {
		tmpprunestep=prunestep(ivariable,X,ivalue,ixindex,ivarindex,idepth,ioutcome)
		tmppruneMH=MHprune(probgrow,probprune,alpha,beta,initsigma,sigma_mu,tmpprunestep$num_term,tmpprunestep$num_pred,tmpprunestep$num_value,tmpprunestep$num_2g_int_nodes,tmpprunestep$num_obs_in_pick_node,tmpprunestep$num_obs_in_left_child,tmpprunestep$num_obs_in_right_child,tmpprunestep$out_in_left_child,tmpprunestep$out_in_right_child,tmpprunestep$out_in_pick_node,tmpprunestep$depth_pick)
		if(tmppruneMH$accept==1){
			tmpvalue[[1]]=tmpprunestep$value
			tmpxindex[[1]]=tmpprunestep$xindex
			tmpdepth[[1]]=tmpprunestep$depth
			tmpvarindex[[1]]=tmpprunestep$varindex
			tmpvariable[[1]]=tmpprunestep$variable	
		}else{
			tmpvalue[[1]]=ivalue
			tmpxindex[[1]]=ixindex
			tmpdepth[[1]]=idepth
			tmpvarindex[[1]]=ivarindex
			tmpvariable[[1]]=ivariable
		}
	} 
	
	tree=tmpvariable[[1]]

	if(is.vector(tree)){
		tmpnode=rnorm(1,mean=(sigma_mu^2*sum(ioutcome))/(sum(tree)*sigma_mu^2+initsigma^2),sd=sqrt((initsigma^2*sigma_mu^2)/(initsigma^2+sum(tree)*sigma_mu^2)))
		tmp=rep(tmpnode,length(tree))	
	} else {
		tmp=numeric(length(transy))
		treesize=ncol(tree)
		tmpnorm=rnorm(treesize)
		tmpnode=numeric(treesize)
		for(k in 1:treesize){
			tmpnode[k]=(sigma_mu^2*sum(ioutcome[tree[,k]==1]))/(sum(tree[,k])*sigma_mu^2+initsigma^2)+tmpnorm[k]*sqrt((initsigma^2*sigma_mu^2)/(initsigma^2+sum(tree[,k])*sigma_mu^2))
			tmp[tree[,k]==1]=tmpnode[k]
		}
	}
	
	mu_ij_list[[1]]=tmpnode
	tmpmu[,1]=tmp

	for(tree_count in 2:m){
		if(tree_count==(m-1)){
			ioutcome=transy-rowSums(tmpmu[,1:(tree_count-1)])-mu_ij[,m]
		} else if (tree_count==m){
			ioutcome=transy-rowSums(tmpmu[,1:(tree_count-1)])
		} else if (tree_count==2){
			ioutcome=transy-tmpmu[,1]-rowSums(mu_ij[,3:m])
		} else {
			ioutcome=transy-rowSums(tmpmu[,1:(tree_count-1)])-rowSums(mu_ij[,(tree_count+1):m])
		}
		
		ivalue=value[[tree_count]]
		ixindex=xindex[[tree_count]]
		idepth=depth[[tree_count]]
		ivarindex=varindex[[tree_count]]
		ivariable=variable[[tree_count]]

		proposal_type=sample(2,1)

		if(proposal_type==1){
			tmpgrowstep=growstep(ivariable,X,ivalue,ixindex,ivarindex,idepth,ioutcome)
			tmpgrowMH=MHgrow(probgrow,probprune,alpha,beta,initsigma,sigma_mu,tmpgrowstep$num_term,tmpgrowstep$num_pred,tmpgrowstep$num_value,tmpgrowstep$num_2g_int_nodes,tmpgrowstep$num_obs_in_pick_node,tmpgrowstep$num_obs_in_left_child,tmpgrowstep$num_obs_in_right_child,tmpgrowstep$out_in_left_child,tmpgrowstep$out_in_right_child,tmpgrowstep$out_in_pick_node,tmpgrowstep$depth_pick)
	
			if(tmpgrowMH$accept==1){
				tmpvalue[[tree_count]]=tmpgrowstep$value
				tmpxindex[[tree_count]]=tmpgrowstep$xindex
				tmpdepth[[tree_count]]=tmpgrowstep$depth
				tmpvarindex[[tree_count]]=tmpgrowstep$varindex
				tmpvariable[[tree_count]]=tmpgrowstep$variable	
			}else{
				tmpvalue[[tree_count]]=ivalue
				tmpxindex[[tree_count]]=ixindex
				tmpdepth[[tree_count]]=idepth
				tmpvarindex[[tree_count]]=ivarindex
				tmpvariable[[tree_count]]=ivariable
			}
		} else {
			tmpprunestep=prunestep(ivariable,X,ivalue,ixindex,ivarindex,idepth,ioutcome)
			tmppruneMH=MHprune(probgrow,probprune,alpha,beta,initsigma,sigma_mu,tmpprunestep$num_term,tmpprunestep$num_pred,tmpprunestep$num_value,tmpprunestep$num_2g_int_nodes,tmpprunestep$num_obs_in_pick_node,tmpprunestep$num_obs_in_left_child,tmpprunestep$num_obs_in_right_child,tmpprunestep$out_in_left_child,tmpprunestep$out_in_right_child,tmpprunestep$out_in_pick_node,tmpprunestep$depth_pick)
			if(tmppruneMH$accept==1){
				tmpvalue[[tree_count]]=tmpprunestep$value
				tmpxindex[[tree_count]]=tmpprunestep$xindex
				tmpdepth[[tree_count]]=tmpprunestep$depth
				tmpvarindex[[tree_count]]=tmpprunestep$varindex
				tmpvariable[[tree_count]]=tmpprunestep$variable	
			}else{
				tmpvalue[[tree_count]]=ivalue
				tmpxindex[[tree_count]]=ixindex
				tmpdepth[[tree_count]]=idepth
				tmpvarindex[[tree_count]]=ivarindex
				tmpvariable[[tree_count]]=ivariable
			}
		} 

		tree=tmpvariable[[tree_count]]

		if(is.vector(tree)){
			tmpnode=rnorm(1,mean=(sigma_mu^2*sum(ioutcome))/(sum(tree)*sigma_mu^2+initsigma^2),sd=sqrt((initsigma^2*sigma_mu^2)/(initsigma^2+sum(tree)*sigma_mu^2)))
			tmp=rep(tmpnode,length(tree))	
		} else {
			tmp=numeric(length(transy))
			treesize=ncol(tree)
			tmpnorm=rnorm(treesize)
			tmpnode=numeric(treesize)
			for(k in 1:treesize){
				tmpnode[k]=(sigma_mu^2*sum(ioutcome[tree[,k]==1]))/(sum(tree[,k])*sigma_mu^2+initsigma^2)+tmpnorm[k]*sqrt((initsigma^2*sigma_mu^2)/(initsigma^2+sum(tree[,k])*sigma_mu^2))
				tmp[tree[,k]==1]=tmpnode[k]
			}
		}
		mu_ij_list[[tree_count]]=tmpnode
		tmpmu[,tree_count]=tmp
		
	}
	
	tmppredy=rowSums(tmpmu)

	return(list(un_trans=(tmppredy-beta0)/beta1,value=tmpvalue,xindex=tmpxindex,depth=tmpdepth,varindex=tmpvarindex,variable=tmpvariable,mu_ij=tmpmu,mu_ij_list=mu_ij_list))

}

bartdraw=function(y,X,base,initsigma,sigma_mu,m,alpha,beta,probgrow,probprune,numpost=1000,numburnin=100){

	beta1=(2*base)/(max(y)-min(y))
	beta0=base-max(y)*beta1
	transy=beta0+beta1*y
	sigma=beta1*initsigma
	init_mu_ij=matrix(1/m*mean(transy),length(transy),m)
	init_value=list()
	init_xindex=list()
	init_depth=list()
	init_varindex=list()
	init_variable=list()
	for(i in 1:m){
		init_value[[i]]=NA
		init_xindex[[i]]=NA
		init_depth[[i]]=1
		init_varindex[[i]]=1
		init_variable[[i]]=rep(1,length(transy))	
	}
	onedraw=bart1draw(transy,init_mu_ij,sigma,sigma_mu,X,init_value,init_xindex,init_depth,init_varindex,init_variable,m,alpha,beta,probgrow,probprune,beta0,beta1)
	for(loop in 2:(numburnin+numpost)) onedraw=bart1draw(transy,onedraw$mu_ij,sigma,sigma_mu,X,onedraw$value,onedraw$xindex,onedraw$depth,onedraw$varindex,onedraw$variable,m,alpha,beta,probgrow,probprune,beta0,beta1)
	
	return(onedraw$un_trans)
}

keep <- function(..., x = c()){
    	if (length(x) > 0){
      		if (!is.character(x)) stop("x must contain character vector")
      		L <- ls(name = parent.frame())
      		rm(list = L[!L %in% x], pos = parent.frame())
      		return(invisible(ls(name = parent.frame())))
    	}
    	dots <- match.call(expand.dots = FALSE)$...
    	if (length(dots) && !all(sapply(dots, function(x) is.symbol(x) || 
                                    is.character(x)))) 
      	stop("... must contain names or character strings")
    
    	names <- sapply(dots, as.character)
    	L <- ls(name = parent.frame())
    	rm(list = L[!L %in% names], pos = parent.frame())
    
    	return(invisible(ls(name = parent.frame())))
    
}