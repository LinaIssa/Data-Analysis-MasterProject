! Wrapper to dgesvx. Solve is 'in-place' (right-hand side in input is solution on output)
! =========================================
SUBROUTINE SQUARE_SOLVE(A,neq,XX)
! =========================================
  implicit none
  integer::neq
  logical::invert
  character(1)                :: FACT,TRANS,EQUED
  integer                     :: nrhs,ldab,ldafb,ldb,ldx,lda
  real,dimension(neq,neq)  :: A
  real,dimension(neq,neq):: AFB
  integer,dimension(neq)       :: IPIV
  real,dimension(neq)       :: Rs,Cs
  real,dimension(neq)::XX,B
  real                        :: Rcond,FERR,BERR
  real,    dimension(4*neq)     :: WORK
  integer, dimension(neq)     :: IWORK
  integer                     :: INFO
  print*,'neq=',neq
  print*,'A(1,:)',A(1,:)
  B=xx
  nrhs=1
  lda=neq
  ldafb=neq
  ldb=neq
  ldx=neq
  FACT = 'E' ! Equilibrate if necessary
  EQUED= 'N' ! Say that matrix has not yet been equilibrated. 
  TRANS= 'N' ! no transpose
  call DGESVX( FACT, TRANS, neq, NRHS, A, LDA, AFB, LDAFB, &
       IPIV, EQUED, Rs, Cs, B, LDB, XX, LDX, RCOND, FERR, BERR, &
       WORK, IWORK, INFO )
  if(info<0) then
     print*,'Solve matrix: Illegal value for element',-info
  endif
  if (info>0) then
     if (info==neq+1) then
        print*,'Warning, matrix singular to working precision..(EQUED=',EQUED,')'
     else
        print*,'Solve matrix, singular element at ',INFO, '(EQUED=',EQUED,')'
     end if
  endif
end SUBROUTINE 
