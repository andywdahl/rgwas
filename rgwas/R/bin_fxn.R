bin_fxn	<- function( mu0, sd0, zeros ){ # O(N)

	z		<- -mu0/sd0

	l.binp	<- l.binp_fxn( z, zeros )

	psi					<- exp( dnorm( z, log=T ) - l.binp )
	psi[-zeros]	<- -psi[-zeros]
	muY					<- mu0 - sd0*psi
	SigY				<- sd0^2*(1 - z*psi - psi^2)

	list( l.binp=l.binp, muY=muY, SigY=SigY )
}

l.binp_fxn	<- function( z, zeros ){
	l.binp	<- rep( NA, length(z) )
	l.binp[-zeros]<- pnorm(	z[-zeros]	, lower.tail=F, log.p=T )
	l.binp[zeros]	<- pnorm(	z[zeros]	, lower.tail=T, log.p=T )
	l.binp
}
