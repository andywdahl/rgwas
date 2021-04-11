score_K_ps <- function( Yb, Yq, G, X=NULL, K, n.folds=10, mc.cores=1, nrun=10, tau=0, ... )
	scores	<- as.numeric(parallel::mclapply( 1:n.folds, mc.cores=mc.cores, function(kk)
{
	cat( '\nrunning fold ', kk, '....\n' )
	sub1	<- sample( nrow(G), nrow(G)/2, replace=F )

	if( !is.null(X) ){
		X1 	<- X [sub1 ,,drop=F]
		X2 	<- X [-sub1,,drop=F]
	} else {
		X1 	<- NULL
		X2 	<- NULL
	}

	fit1	<- mfmr( Yb[ sub1,,drop=F], Yq[ sub1,,drop=F], G[ sub1,,drop=F], X1, K, seeds=1:nrun+898+737*kk				, nrun=nrun, ... )
	fit2	<- mfmr( Yb[-sub1,,drop=F], Yq[-sub1,,drop=F], G[-sub1,,drop=F], X2, K, seeds=1:nrun+898+737*kk + 89938, nrun=nrun, ... )

	# original clusters
	#z1		<- apply( fit1$pmat		, 1, which.max )
	#z2		<- apply( fit2$pmat		, 1, which.max )
	z1		<- apply( fit1$pmat		, 1, function(x) ifelse( max(x)>tau, which.max(x), NA ) )
	z2		<- apply( fit2$pmat		, 1, function(x) ifelse( max(x)>tau, which.max(x), NA ) )

	# new clusters
	newpmat1	<- fit_z( Yb[ sub1,,drop=F], Yq[ sub1,,drop=F], G[ sub1,,drop=F], X1, fit2$out )
	newpmat2	<- fit_z( Yb[-sub1,,drop=F], Yq[-sub1,,drop=F], G[-sub1,,drop=F], X2, fit1$out )
	rm( fit1, fit2 )

	# convert to matrix 
	newz1	<- apply( newpmat1, 1, function(x) ifelse( max(x)>tau, which.max(x), NA ) )
	newz2	<- apply( newpmat2, 1, function(x) ifelse( max(x)>tau, which.max(x), NA ) )
	Z1	<- model.matrix( ~ -1 + as.factor(newz1) )
	Z2	<- model.matrix( ~ -1 + as.factor(newz2) )
	#Z1	<- newpmat1
	#Z2	<- newpmat2

	# compare new/original clusters
	1/2*(ps_dist( Z1, z1 )+ps_dist( Z2, z2 ))
}))

ps_dist	<- function( Z, z0 )
	min(sapply( 1:ncol(Z), function(k)
{
	Zk	<- Z[which( z0==k ),]
	nk	<- nrow(Zk)
	#D1k	<- Zk %*% t(Zk)
	#1/nk/(nk-1) * ( sum( D1k )/2 - sum(diag( D1k )) )
	( sum( colSums( Zk )^2 ) - nk )/(nk^2-nk)
}))
