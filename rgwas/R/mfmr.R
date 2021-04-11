mfmr	<- function(
	Yb, Yq, G, X=NULL, K,
	nrun=10, init_sd=1e-2, seeds=1:nrun+898,
	mc.cores=1,
	trace=FALSE,
	...
){

	if( !is.null(Yb) ){
		stopifnot( all( Yb %in% 0:1 ) )
		stopifnot( length( which( Yb[,1] == 0 ) ) > 3 )
	}
	stopifnot( all( !is.na( c(Yb,Yq,G,X) ) ) )
	
	if( var(G[,1]) > 1e-8 )
		warning( 'G should almost always include an intercept column' )

	outs	<- parallel::mclapply( 1:nrun, mc.cores=mc.cores, function(ii){
		if( nrun > 1 ) cat( 'Running random run #', ii )
		ipar	<- initialize_fxn( Yb, Yq, G, X, K, init_sd, seed=seeds[ii], ... )
		out	<- list( ll=NA, conv=F )
		tryCatch({
			out	<- mfmr_em_alg( Yb, Yq, G, X, K, ifelse( ii == 1 | mc.cores == 1, trace, FALSE ), ipar$N, ipar$Q, ipar$S, ipar$P, ipar$B, ipar$p, ipar$al, ipar$be, ipar$Si, ... )
			if( nrun > 1 ) cat( ', ll=', out$ll, '\n' )
		}, error=function(e)
			print(e)
		)
			#cat( 'Error in an MFMR run, ID=', ii, '; Error: ', e, '\n' )
		out
	})

	lls		<- sapply( outs, function(out) out$ll )
	if( any(is.na(lls)) )
		print( paste( 'Warning: ', sum(is.na(lls)), 'runs failed' ) )
	out		<- outs[[which.min(lls)]]
	rm( outs )

	pmat	<- fit_z( Yb, Yq, G, X, out )

	list( pmat=pmat, out=out, ll=max(lls), all_lls=lls )
}
