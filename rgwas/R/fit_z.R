fit_z	<- function( Yb, Yq, G, X, out ){

	if( is.null( Yb ) )
		return( fit_z_nobin( Yq, G, X, out ) )

	N	<- nrow(Yb)
	B	<- ncol(Yb)
	K	<- length(out$pvec)

	if( is.null( X ) ){
		Xalpha	<- matrix( 0, N, B+ncol(Yq) )
	} else {
		Xalpha	<- X %*% out$alpha
	}

	beta	<- out$beta
	Lambda_qq	<- solve( out$Sigma[-(1:B),-(1:B)] )

	pmat	<- matrix( NA, N, K )
	for( k in 1:K ){

		qerr			<- Yq - Xalpha[,-(1:B)] - G %*% beta[k,,-(1:B)]
		Lam_qerr	<- qerr %*% Lambda_qq
		mu0				<- Xalpha[,1:B]	+ G %*% beta[k,,1:B] + Lam_qerr %*% out$Sigma[-(1:B),1:B]

		l.binp	<- matrix( NA, N, B )
		for( b in 1:B )
			l.binp[,b]	<- bin_fxn( mu0[,b], sd0=1, zeros=which(Yb[,b]==0) )$l.binp
		pmat[,k]	<- log(out$pvec[k]) + rowSums(l.binp) - 1/2*rowSums( Lam_qerr * qerr )
	}
	lpmat2pmat( pmat )
}

fit_z_nobin	<- function( Y, G, X, out ){

	N	<- nrow(Y)
	K	<- length(out$pvec)

	if( is.null( X ) ){
		Xalpha	<- matrix( 0, N, ncol(Y) )
	} else {
		Xalpha	<- X %*% out$alpha
	}

	Lambda	<- solve( out$Sigma )

	pmat	<- matrix( NA, N, K )
	for( k in 1:K ){
		qerr			<- Y - Xalpha - G %*% out$beta[k,,]
		Lam_qerr	<- qerr %*% Lambda
		pmat[,k]	<- log(out$pvec[k]) - 1/2*rowSums( Lam_qerr * qerr )
	}
	lpmat2pmat( pmat )
}
