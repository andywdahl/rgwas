mfmr_em_alg_noquant	<- function(
	Yb, G, X, K, trace,
	N, Q, S, B,
	pvec, alpha, beta,
	maxit=1e3, tol=rep(1e-3,2)
){

	al_eps	<- 0
	Xalpha	<- matrix( 0, N, B )
	if( Q > 0 ){
		Xalpha	<- X %*% alpha
		XtXiXt	<- solve(t(X)%*%X) %*% t(X)
	}

	if( trace )
		cat(sprintf("%10d, %10.3e, %10.3e, %10.3e, %10.3e, %10.3e, %4f, %4f", 0, NA, NA, NA, beta[1,1,1], beta[2,1,1], pvec[1], pvec[2] ),'\n')
	for( it in 1:maxit ){
		olda		<- alpha
		oldb		<- beta

		pmat	<- matrix( NA, N, K )
		muY		<- array( NA, dim=c(N,K,B) )
		for( k in 1:K ){

			mu0	<- Xalpha[,1:B]	+ G %*% beta[k,,1:B]

			l.binp	<- matrix( NA, N, B )
			for( b in 1:B ){
				binout			<- bin_fxn( mu0[,b], sd0=1, zeros=which(Yb[,b]==0) )
				l.binp[,b]	<- binout$l.binp
				muY	[,k,b]	<- binout$muY
				rm( binout )
			}
			pmat[,k]	<- log(pvec[k]) + rowSums(l.binp)
		}
 		pmat	<- lpmat2pmat( pmat )
		pvec	<- colMeans( pmat )

		if( !is.null(X) ){
			meanmat	<- apply( muY, 3, function(muY_b) rowSums( muY_b * pmat ) )
			delta		<- meanmat - sapply( 1:B, function(p) rowSums( G * (pmat %*% beta[,,p]) ))
			alpha		<- XtXiXt %*% delta
			Xalpha	<- X %*% alpha
			al_eps	<- sqrt( mean((alpha-olda)^2)/mean(alpha^2) )
			if( Q == 1 ) alpha		<- matrix( alpha,		Q, B )
		}

		for( k in 1:K ){
			xG				<- apply( G, 2, function(g) g*pmat[,k] )
			beta[k,,]	<- solve( t(G)%*%xG ) %*% t(xG)%*%as.matrix( muY[,k,] - Xalpha )
		}
		be_eps	<- sqrt( mean((beta-oldb)^2)/ mean( sapply( 1:K, function(k) c( beta[k,,] - apply( beta, 2:3, mean ) )^2 ) ) )

		errmat	<- 0
		for( k in 1:K ){
			sqrt.err.k	<- apply( muY[,k,] - Xalpha - G%*%beta[k,,], 2, function(e) e*sqrt(pmat[,k]) )
			errmat			<- errmat + t(sqrt.err.k) %*% sqrt.err.k
		}

		if( trace )
			cat(sprintf("%10d, %10.3e, %10.3e, %10.3e, %10.3e, %4f, %4f", it, al_eps, be_eps, beta[1,1,1], beta[2,1,1], pvec[1], pvec[2] ),'\n')
		if( it > 10 & all( c( al_eps, be_eps ) < tol, na.rm=T ) ) break

	}
	ll	<- ll_fxn_noquant( Yb, Xalpha, G, beta, pvec, N, B, K )
	list( pvec=pvec, alpha=alpha, beta=beta, Sigma=NA, ll=ll, niter=it, conv=it<maxit )
}

ll_fxn_noquant	<- function( Yb, Xalpha, G, beta, pvec, N, B, K ){

	all_lls_NxK	<- sapply( 1:K, function(k){
		mu0			<- Xalpha[,1:B]	+ G %*% beta[k,,1:B]
		log(pvec[k]) + rowSums( sapply( 1:B, function(b) l.binp_fxn( z = - mu0[,b], zeros=which(Yb[,b]==0) ) ))
	})
	offsets		<- apply( all_lls_NxK, 1, max )
	all_lls_N	<- sum(log(rowSums(exp( all_lls_NxK - (offsets %o% rep(1,K)) )))) + sum(offsets)

	-N*B/2*log(2*pi) + all_lls_N
}
