droptest	<- function(
	Yb, Yq, G, X=NULL, test_inds, K,
	mc.cores=1,
	mfmrx=FALSE, pmat=NULL,
	...
){

	if( !is.null(Yb) & is.null( colnames(Yb) ) ) 	colnames(Yb)	<- paste0( 'Yb'	, 1:ncol(Yb))
	if( !is.null(Yq) & is.null( colnames(Yq) ) ) 	colnames(Yq)	<- paste0( 'Yq'	, 1:ncol(Yq))
	if( !is.null(G ) & is.null( colnames(G ) ) )	colnames(G)	<- paste0( 'G'	, 1:ncol(G) )

	phens	<- c( colnames( Yb ), colnames( Yq ) )
	P	<- length(phens)

	##### set up fixed out
	if( mfmrx & is.null(pmat) ){
		mfmr_time	<- system.time({
			pmat	<- mfmr( Yb, Yq, G, X, K, mc.cores=mc.cores, ... )$pmat
		})[3]
		cat( 'mfmr time:', round( mfmr_time/60, 1 ), 'min.', '\n' )
	}

	##### test each column of G
	gxe_outs	<- parallel::mclapply( 1:ncol(G), mc.cores=mc.cores, function(test_ind){
		if( ! test_ind %in% test_inds ) return(NA)
		cat( 'Running G index ', test_ind, '\n' )

		##### fit mfmr
		mc.cores.mfmr	<- ifelse( mc.cores <= length( test_inds ), 1, floor(mc.cores/length(test_inds)) )
		if( ! mfmrx  ){
			mfmr_time	<- system.time({
				pmat	<- mfmr( Yb, Yq, G[,-test_ind,drop=F], cbind( X, G[,test_ind] ), K, mc.cores=mc.cores.mfmr, ... )$pmat
			})[3]
			cat( 'mfmr time:', round( mfmr_time/60, 1 ), 'min.', '\n' )
		}


		##### test
		Gxpmat	<- t(sapply( 1:nrow(G), function(i) G[i,-test_ind,drop=F] %x% pmat[i,] ))
		out	<- lapply( 1:P, function(pp){
			bin	<- ifelse( is.null(Yb), FALSE, (pp <= ncol(Yb)) )
			y		<- cbind( Yb, Yq )[,pp]
			out1	<- interxn_test( X=cbind( X, pmat[,-K,drop=F], Gxpmat ), y=y, g=G[,test_ind], pmat=pmat[,-K,drop=F], bin=bin )
			list( pvals=out1$pvals, summ=summary( out1$fit )$coef )
		})
		out
	})

	betas	<- array( NA, dim=c(K,ncol(G),P), dimnames=list( 1:K											, colnames(G), phens ) )
	ses		<- array( NA, dim=c(K,ncol(G),P), dimnames=list( 1:K											, colnames(G), phens ) )
	pvals	<- array( NA, dim=c(3,ncol(G),P), dimnames=list( c( 'Hom', 'Het', 'Glob' ), colnames(G), phens ) )
	for( s in test_inds )
		for( p in 1:P )
	{
		out_sp	<- gxe_outs[[s]][[p]]
		betas	[,s,p]	<- out_sp$summ[paste0('gxe',1:K),'Estimate']
		ses		[,s,p]	<- out_sp$summ[paste0('gxe',1:K),'Std. Error']
		pvals	[,s,p]	<- out_sp$pvals
	}
	return(list( betas=betas, ses=ses, pvals=pvals ))
}
