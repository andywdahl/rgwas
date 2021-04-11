score_K	<- function( Yb, Yq, G, X=NULL, K, n.folds=10, mc.cores=1, nrun=10, ... ){

	B		<- ncol(Yb)
	P		<- ncol(cbind(Yb,Yq))
	N		<- nrow(cbind(Yb,Yq))
	subs<- split( 1:N, sample( rep(1:n.folds,N), N ) )

	score_matrix	<- matrix( NA, n.folds, 2 )
	colnames(score_matrix)	<- c( 'oos_ll', 'negentropy' )
	cat( 'running fold ' )
	for( kk in 1:n.folds ){
		cat( kk, '...' )
		sub			<- subs[[kk]]
		
		trn_Yb	<- Yb[-sub,,drop=F]
		trn_Yq	<- Yq[-sub,,drop=F]
		trn_G 	<- G [-sub,,drop=F]

		tst_Yb	<- Yb[sub,,drop=F]
		tst_Yq	<- Yq[sub,,drop=F]
		tst_G 	<- G [sub,,drop=F]
		Ntest	<- nrow(cbind(tst_Yb,tst_Yq))

		if( !is.null(X) ){
			trn_X 	<- X [-sub,,drop=F]
			tst_X 	<- X [sub,,drop=F]
		} else {
			trn_X 	<- NULL
			tst_X 	<- NULL
		}

		if( K > 1 ){
			out	<- mfmr( trn_Yb, trn_Yq, trn_G, trn_X, K, mc.cores=mc.cores, nrun=nrun, ... )$out
			if( !is.null(X) ){
				Xalpha	<- tst_X %*% out$alpha
			} else {
				Xalpha 	<- matrix( 0, Ntest, P )
			}
			score_matrix[kk,1]	<- ll_fxn(		tst_Yb, tst_Yq, Xalpha, tst_G, out$beta	, out$Sigma, out$pvec	, Ntest, P, K, B )
			score_matrix[kk,2]	<- -sum( out$pvec*log(out$pvec) )
		} else {
			out	<- mfmr_em_alg_hom( trn_Yb, trn_Yq					, X=cbind( trn_G, trn_X ), nrun=nrun, ... )
			score_matrix[kk,1]	<- ll_fxn_hom(tst_Yb, tst_Yq, Xalpha=cbind( tst_G, tst_X ) %*% out$alpha, out$Sigma, Ntest, P, B )
		}
	}
	cat( 'Done!\n' )
	score_matrix
}
