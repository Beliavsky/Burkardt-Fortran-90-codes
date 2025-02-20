subroutine amino_print ( acid_num, acid_sym )

!*****************************************************************************80
!
!! amino_print() prints the amino acid parameters.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer ACID_NUM, the number of amino acids.
!
!    character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
  implicit none

  integer acid_num

  integer acid_i
  character ( len = 27 ) acid_name
  character acid_sym(acid_num)
  character c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  Amino Acid Symbol'
  write ( *, '(a)' ) ' '
  do acid_i = 1, acid_num
    c = acid_sym(acid_i)
    call ch_to_amino_name ( c, acid_name )
    write ( *, '(i3,2x,a,2x,a)' ) acid_i, acid_sym(acid_i), acid_name
  end do

  return
end
subroutine binomial_sample ( a, b, x )

!*****************************************************************************80
!
!! binomial_sample() samples the Binomial PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Kennedy, James Gentle,
!    Algorithm BU,
!    Statistical Computing,
!    Dekker, 1980.
!
!  Input:
!
!    integer A, the number of trials.
!    1 <= A.
!
!    real ( kind = rk ) B, the probability of success on one trial.
!    0.0D+00 <= B <= 1.0.
!
!  Output:
!
!    integer X, a sample of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) u
  integer x

  x = 0

  do i = 1, a

    call random_number ( harvest = u )

    if ( u <= b ) then
      x = x + 1
    end if

  end do

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! ch_cap() capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character C, the character to capitalize.
!
!  Output:
!
!    character C, the capitalized character.
!
  implicit none

  character c
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! ch_eqi() is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character C1, C2, the characters to compare.
!
!  Output:
!
!    logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_next ( line, cval, done )

!*****************************************************************************80
!
!! ch_next() "reads" space-separated characters from a string, one at a time.
!
!  Example:
!
!    Input:
!
!      LINE = ' A  B, C    DE  F'
!
!    Output:
!
!      'A', 'B', 'C', 'D', 'E', 'F', and then blanks.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) LINE, a string, presumably containing
!    characters, possibly separated by spaces or commas.
!
!    logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!
!  Output:
!
!    character CVAL.  If DONE is FALSE, then CVAL contains the
!    "next" character read from LINE.  If DONE is TRUE, then
!    CVAL is blank.
!
!    logical DONE: the routine sets DONE to FALSE if another character
!    was read, or TRUE if no more characters could be read.
!
  implicit none

  character cval
  logical done
  integer i
  character ( len = * ) line
  integer, save :: next = 1

  if ( done ) then
    next = 1
    done = .false.
  end if

  do i = next, len ( line )

    if ( line(i:i) /= ' ' .and. line(i:i) /= ',' ) then
      cval = line(i:i)
      next = i + 1
      return
    end if

  end do

  done = .true.
  next = 1
  cval = ' '

  return
end
subroutine ch_to_amino_name ( c, amino_name )

!*****************************************************************************80
!
!! ch_to_amino_name() converts a character to an amino acid name.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl Branden, John Tooze,
!    Introduction to Protein Structure,
!    Garland Publishing, 1991.
!
!  Input:
!
!    character C, the one letter code for an amino acid.
!    Lower and upper case letters are treated the same.
!
!  Output:
!
!    character ( len = * ) AMINO_NAME, the full name of the
!    corresponding amino acid.  The longest name is 27 characters.  
!    If the input code is not recognized, then AMINO_NAME will be 
!    set to '???'.
!
  implicit none

  integer, parameter :: n = 23

  character ( len = * ) amino_name
  character ( len = 27 ), dimension ( n ) :: amino_table = (/ &
    'Alanine                    ', &
    'Aspartic acid or Asparagine', &
    'Cysteine                   ', &
    'Aspartic acid              ', &
    'Glutamic acid              ', &
    'Phenylalanine              ', &
    'Glycine                    ', &
    'Histidine                  ', &
    'Isoleucine                 ', &
    'Lysine                     ', &
    'Leucine                    ', &
    'Methionine                 ', &
    'Asparagine                 ', &
    'Proline                    ', &
    'Glutamine                  ', &
    'Arginine                   ', &
    'Serine                     ', &
    'Threonine                  ', &
    'Valine                     ', &
    'Tryptophan                 ', &
    'Undetermined amino acid    ', &
    'Tyrosine                   ', &
    'Glutamic acid or Glutamine ' /)
  character c
  logical ch_eqi
  character, dimension ( n ) :: c_table = (/ &
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', &
    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', &
    'X', 'Y', 'Z' /)
  integer i

  do i = 1, n
    if ( ch_eqi ( c, c_table(i) ) ) then
      amino_name = amino_table(i)
      return
    end if
  end do

  amino_name = '???'

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! ch_to_digit() returns the value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!  Output:
!
!    integer DIGIT, the corresponding value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine comp_param_print ( acid_num, acid_sym, comp_max, comp_num, beta, &
  beta_sum, comp_weight )

!*****************************************************************************80
!
!! comp_param_print() prints the parameters for the mixture components.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer ACID_NUM, the number of amino acids.
!
!    character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    integer COMP_MAX, the maximum number of Dirichlet 
!    mixture components.
!
!    integer COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
!    real ( kind = rk ) BETA(ACID_NUM,COMP_MAX); BETA(I,J) is the 
!    parameter for the J-th acid in the I-th Dirichlet mixture component.
!
!    real ( kind = rk ) BETA_SUM(COMP_MAX), the sum of the values of
!    BETA(ACID_I,COMP_I) for a given component COMP_I.
!
!    real ( kind = rk ) COMP_WEIGHT(COMP_NUM), the mixture weight of each
!    component.  These values should be nonnegative, and sum to 1.  They
!    represent the relative proportion of each component in the mixture.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer acid_num
  integer comp_max

  integer acid_i
  character acid_sym(acid_num)
  integer comp_i
  real ( kind = rk ) beta(acid_num,comp_max)
  real ( kind = rk ) beta_sum(comp_max)
  integer comp_num
  real ( kind = rk ) comp_weight(comp_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of components = ', comp_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(''Compon:'',20i8)' ) ( comp_i, comp_i = 1, comp_num )
  write ( *, '(''Weight:'',20f8.4)' ) comp_weight(1:comp_num)
  write ( *, '(a)' ) ' '

  do acid_i = 1, acid_num
    write ( *, '(i2,2x,a1,2x,20f8.4)' ) acid_i, acid_sym(acid_i), &
      beta(acid_i,1:comp_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a3,4x,20f8.4)' ) 'Sum', beta_sum(1:comp_num)

  return
end
subroutine dirichlet_mean ( n, a, mean )

!*****************************************************************************80
!
!! dirichlet_mean() returns the means of the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components.
!
!    real ( kind = rk ) A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be positive.
!
!  Output:
!
!    real ( kind = rk ) MEAN(N), the means of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) mean(n)

  if ( any ( a(1:n) < 0.0D00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_MEAN - Fatal error!'
    write ( *, '(a)' ) '  At least one entry of A is negative!'
    stop
  end if

  if ( all ( a(1:n) == 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_MEAN - Fatal error!'
    write ( *, '(a)' ) '  All entries of A are zero!'
    stop
  end if

  mean(1:n) = a(1:n)

  call r8vec_unit_sum ( n, mean )

  return
end
subroutine dirichlet_mix_check ( comp_num, elem_max, elem_num, a, comp_weight )

!*****************************************************************************80
!
!! dirichlet_mix_check() checks the parameters of a Dirichlet mixture PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer COMP_NUM, the number of components in the 
!    Dirichlet mixture density, that is, the number of distinct Dirichlet PDF's
!    that are mixed together.
!
!    integer ELEM_MAX, the leading dimension of A, which
!    must be at least ELEM_NUM.
!
!    integer ELEM_NUM, the number of elements of an 
!    observation.
!
!    real ( kind = rk ) A(ELEM_MAX,COMP_NUM), the probabilities for 
!    element ELEM_NUM in component COMP_NUM.
!    Each A(I,J) should be greater than or equal to 0.0.
!
!    integer COMP_WEIGHT(COMP_NUM), the mixture weights of
!    the densities. These do not need to be normalized.  The weight of a given 
!    component is the relative probability that that component will be used 
!    to generate the sample.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer comp_num
  integer elem_max
  integer elem_num

  real ( kind = rk ) a(elem_max,comp_num)
  integer comp_i
  real ( kind = rk ) comp_weight(comp_num)
  integer elem_i
  logical positive

  do comp_i = 1, comp_num

    do elem_i = 1, elem_num
      if ( a(elem_i,comp_i) < 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIRICHLET_MIX_CHECK - Fatal error!'
        write ( *, '(a)' ) '  A(ELEM,COMP) < 0.'
        write ( *, '(a,i6)' ) '  COMP = ', comp_i
        write ( *, '(a,i6)' ) '  ELEM = ', elem_i
        write ( *, '(a,g14.6)' ) '  A(COMP,ELEM) = ', a(elem_i,comp_i)
        stop
      end if
    end do

  end do

  positive = .false.

  do comp_i = 1, comp_num

    if ( comp_weight(comp_i) < 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIRICHLET_MIX_CHECK - Fatal error!'
      write ( *, '(a)' ) '  COMP_WEIGHT(COMP) < 0.'
      write ( *, '(a,i6)' ) '  COMP = ', comp_i
      write ( *, '(a,g14.6)' ) '  COMP_WEIGHT(COMP) = ', comp_weight(comp_i)
      stop
    else if ( 0.0D+00 < comp_weight(comp_i)  ) then
      positive = .true.
    end if

  end do

  if ( .not. positive ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRICHLET_MIX_CHECK - Fatal error!'
    write ( *, '(a)' ) '  All component weights are zero.'
    stop
  end if

  return
end
subroutine dirichlet_multinomial_pdf ( x, a, b, c, pdf )

!*****************************************************************************80
!
!! dirichlet_multinomial_pdf() evaluates a Dirichlet Multinomial PDF.
!
!  Discussion:
!
!    PDF(X)(A,B,C) = Comb(A,B,X) * ( Gamma(C_Sum) / Gamma(C_Sum+A) )
!      Product ( 1 <= I <= B ) Gamma(C(I)+X(I)) / Gamma(C(I))
!
!    where:
!
!      Comb(A,B,X) is the multinomial coefficient C( A; X(1), X(2), ..., X(B) ),
!      C_Sum = Sum ( 1 <= I <= B ) C(I)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kenneth Lange,
!    Mathematical and Statistical Methods for Genetic Analysis,
!    Springer, 1997, page 45.
!
!  Input:
!
!    integer X(B); X(I) counts the number of occurrences of
!    outcome I, out of the total of A trials.
!
!    integer A, the total number of trials.
!
!    integer B, the number of different possible outcomes
!    on one trial.
!
!    integer C(B); C(I) is the Dirichlet parameter 
!    associated with outcome I.
!
!  Output:
!
!    real ( kind = rk ) PDF, the value of the Dirichlet multinomial PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer b

  integer a
  real ( kind = rk ) c(b)
  real ( kind = rk ) c_sum
  integer i
  real ( kind = rk ) pdf
  real ( kind = rk ) pdf_log
  real ( kind = rk ) r8_gamma_log
  integer x(b)

  c_sum = sum ( c(1:b) )

  pdf_log = - r8_gamma_log ( c_sum + real ( a, kind = rk ) ) &
    + r8_gamma_log ( c_sum ) &
    + r8_gamma_log ( real ( a + 1, kind = rk ) )

  do i = 1, b

    pdf_log = pdf_log + r8_gamma_log ( c(i) + real ( x(i), kind = rk ) ) &
      - r8_gamma_log ( c(i) ) - r8_gamma_log ( real ( x(i) + 1, kind = rk ) )
  end do

  pdf = exp ( pdf_log )

  return
end
subroutine dirichlet_sample ( n, a, x )

!*****************************************************************************80
!
!! dirichlet_sample() samples the Dirichlet PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 169.
!
!  Input:
!
!    integer N, the number of components.
!
!    real ( kind = rk ) A(N), the probabilities for each component.
!    Each A(I) should be nonnegative, and at least one should be
!    positive.
!
!  Output:
!
!    real ( kind = rk ) X(N), a sample of the PDF.  The entries 
!    of X should sum to 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a2
  real ( kind = rk ) b2
  real ( kind = rk ) c2
  integer i

  real ( kind = rk ) x(n)

  a2 = 0.0D+00
  b2 = 1.0D+00

  do i = 1, n
    c2 = a(i)
    call gamma_sample ( a2, b2, c2, x(i) )
  end do
!
!  Rescale the vector to have unit sum.
!
  call r8vec_unit_sum ( n, x )

  return
end
subroutine discrete_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! discrete_cdf_inv() inverts the Discrete CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    integer A, the number of probabilities assigned.
!
!    real ( kind = rk ) B(A), the relative probabilities of outcomes 
!    1 through A.  Each entry must be nonnegative.
!
!  Output:
!
!    integer X, the corresponding argument for which
!    CDF(X-1) < CDF <= CDF(X)
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a

  real ( kind = rk ) b(a)
  real ( kind = rk ) b_sum
  real ( kind = rk ) cdf
  real ( kind = rk ) cum
  integer j
  integer x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DISCRETE_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  b_sum = sum ( b(1:a) )

  cum = 0.0D+00

  do j = 1, a

    cum = cum + b(j) / b_sum

    if ( cdf <= cum ) then
      x = j
      return
    end if

  end do

  x = a

  return
end
subroutine discrete_sample ( a, b, x )

!*****************************************************************************80
!
!! discrete_sample() samples the Discrete PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer A, the number of probabilities assigned.
!
!    real ( kind = rk ) B(A), the relative probabilities of 
!    outcomes 1 through A.  Each entry must be nonnegative.
!
!  Output:
!
!    integer X, a sample of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a

  real ( kind = rk ) b(a)
  real ( kind = rk ) cdf
  integer x

  call random_number ( harvest = cdf )

  call discrete_cdf_inv ( cdf, a, b, x )

  return
end
subroutine event_process ( acid_num, alpha, beta, comp_num, p, p_hat, &
  site_num, x_sample )

!*****************************************************************************80
!
!! event_process() updates the mixture weight distribution parameters.
!
!  Discussion:
!
!    This routine updates the values of ALPHA.  It does this by
!    considering the results of the most recent event.  If we knew
!    which component PDF had generated the event, then we would 
!    simply add 1 to the ALPHA for that component.  Instead, we
!    use Bayesian analysis to estimate the proportion of the event
!    that is to be added to each ALPHA.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    BS Everitt, DJ Hand,
!    Finite Mixture Distributions,
!    Chapman and Hall, 1981.
!
!    AFM Smith, UE Makov,
!    A Quasi-Bayes Sequential Procedure for Mixtures,
!    Journal of the Royal Statistical Society,
!    Volume 40, Number 1, B, 1978, pages 106-112.
!
!  Input:
!
!    integer ACID_NUM, the number of amino acids.
!
!    real ( kind = rk ) ALPHA(COMP_NUM), the Dirichlet parameters
!    for the weights.
!
!    real ( kind = rk ) BETA(ACID_NUM,COMP_MAX); BETA(I,J) is the 
!    multinomial Dirichlet parameter for the J-th acid in the I-th Dirichlet 
!    mixture component.
!
!    integer COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
!    real ( kind = rk ) P(COMP_NUM); P(I) is the Bayesian 
!    posterior probability of component I, given the observation of the most 
!    recent event, which is proportional to the probability of the event under 
!    the component I PDF, times the prior probability of component I.
!
!    real ( kind = rk ) P_HAT(COMP_NUM), the prior probablities 
!    of the components.
!
!    integer SITE_NUM, the number of sites observed for 
!    this event.  This value might change from call to call, although in the 
!    demonstration it is kept fixed.
!
!    integer X_SAMPLE(ACID_NUM), the "current event", 
!    namely, the count vector for the number of occurrences of each acid out 
!    of the total of SITE_NUM sites analyzed.  This is the evidence used to 
!    update the "theory" for the value of ALPHA.
!
!  Output:
!
!    real ( kind = rk ) ALPHA(COMP_NUM), the Dirichlet parameters
!    for the weights.
!
!    real ( kind = rk ) P(COMP_NUM); P(I) is the Bayesian 
!    posterior probability of component I, given the observation of the most 
!    recent event, which is proportional to the probability of the event under 
!    the component I PDF, times the prior probability of component I.
!
!    real ( kind = rk ) P_HAT(COMP_NUM), the prior probablities of the components.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer acid_num
  integer comp_num

  real ( kind = rk ) alpha(comp_num)
  real ( kind = rk ) alpha_sum
  real ( kind = rk ) beta(acid_num,comp_num)
  integer comp_i
  real ( kind = rk ) comp_pdf
  real ( kind = rk ) p(comp_num)
  real ( kind = rk ) p_hat(comp_num)
  integer site_num
  integer x_sample(acid_num)
!
!  Sum the parameters.
!
  alpha_sum = sum ( alpha(1:comp_num) )
!
!  Update P_HAT.
!
  do comp_i = 1, comp_num
    p_hat(comp_i) = ( ( alpha_sum - 1.0D+00 ) * p_hat(comp_i) + p(comp_i) ) &
      / alpha_sum
  end do
!
!  Generate the new P's.
!  P(COMP_I) = the Bayesian posterior probability of component I,
!  given the observation of event EVENT_I, which is proportional
!  to the probability of event EVENT_I in the component I PDF,
!  times the prior probability of component I.
!
  do comp_i = 1, comp_num

    call dirichlet_multinomial_pdf ( x_sample, site_num, acid_num, &
      beta(1,comp_i), comp_pdf )

    p(comp_i) = comp_pdf * p_hat(comp_i)

  end do
!
!  Normalize the P's.
!
  if ( sum ( p(1:comp_num) ) == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EVENT_PROCESS - Fatal error!'
    write ( *, '(a)' ) '  The P''s sum to 0.'
    stop
  end if

  call r8vec_unit_sum ( comp_num, p )
!
!  Update the alpha's by adding adding appropriate portions of
!  the most recent event to each component's parameter.
!
  alpha(1:comp_num) = alpha(1:comp_num) + p(1:comp_num)

  return
end
subroutine exponential_01_cdf_inv ( cdf, x )

!*****************************************************************************80
!
!! exponential_01_cdf_inv() inverts the Exponential 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!  Output:
!
!    real ( kind = rk ) X, the corresponding argument.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cdf
  real ( kind = rk ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXPONENTIAL_01_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = - log ( 1.0D+00 - cdf )

  return
end
subroutine exponential_01_sample ( x )

!*****************************************************************************80
!
!! exponential_01_sample() samples the Exponential PDF with parameter 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) X, a sample of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cdf
  real ( kind = rk ) x

  call random_number ( harvest = cdf )

  x = - log ( 1.0D+00 - cdf )

  return
end
subroutine exponential_cdf_inv ( cdf, a, b, x )

!*****************************************************************************80
!
!! exponential_cdf_inv() inverts the Exponential CDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) CDF, the value of the CDF.
!    0.0 <= CDF <= 1.0.
!
!    real ( kind = rk ) A, B, the parameters of the PDF.
!    0.0 < B.
!
!  Output:
!
!    real ( kind = rk ) X, the corresponding argument.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) cdf
  real ( kind = rk ) x

  if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EXPONENTIAL_CDF_INV - Fatal error!'
    write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
    stop
  end if

  x = a - b * log ( 1.0D+00 - cdf )

  return
end
subroutine exponential_sample ( a, b, x )

!*****************************************************************************80
!
!! exponential_sample() samples the Exponential PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, the parameters of the PDF.
!    0.0D+00 < B.
!
!  Output:
!
!    real ( kind = rk ) X, a sample of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) cdf
  real ( kind = rk ) x

  call random_number ( harvest = cdf )

  call exponential_cdf_inv ( cdf, a, b, x )

  return
end
subroutine gamma_sample ( a, b, c, x )

!*****************************************************************************80
!
!! gamma_sample() samples the Gamma PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    Original FORTRAN77 version by Joachim Ahrens, Ulrich Dieter.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Joachim Ahrens, Ulrich Dieter,
!    Generating Gamma Variates by a Modified Rejection Technique,
!    Communications of the ACM, 
!    Volume 25, Number 1, January 1982, pages 47 - 54.
!
!    Joachim Ahrens, Ulrich Dieter,
!    Computer Methods for Sampling from Gamma, Beta, Poisson and
!    Binomial Distributions.
!    Computing, 
!    Volume 12, 1974, pages 223 - 246.
!
!    Joachim Ahrens, KD Kohrt, Ulrich Dieter,
!    Algorithm 599,
!    ACM Transactions on Mathematical Software,
!    Volume 9, Number 2, June 1983, pages 255-257.
!
!  Input:
!
!    real ( kind = rk ) A, B, C, the parameters of the PDF.
!    0.0D+00 < B, 
!    0.0D+00 < C.
!
!  Output:
!
!    real ( kind = rk ) X, a sample of the PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ), parameter :: a1 =   0.3333333D+00
  real ( kind = rk ), parameter :: a2 = - 0.2500030D+00
  real ( kind = rk ), parameter :: a3 =   0.2000062D+00
  real ( kind = rk ), parameter :: a4 = - 0.1662921D+00
  real ( kind = rk ), parameter :: a5 =   0.1423657D+00
  real ( kind = rk ), parameter :: a6 = - 0.1367177D+00
  real ( kind = rk ), parameter :: a7 =   0.1233795D+00
  real ( kind = rk ) b
  real ( kind = rk ) bcoef
  real ( kind = rk ) c
  real ( kind = rk ) co
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ), parameter :: e1 = 1.0D+00
  real ( kind = rk ), parameter :: e2 = 0.4999897D+00
  real ( kind = rk ), parameter :: e3 = 0.1668290D+00
  real ( kind = rk ), parameter :: e4 = 0.0407753D+00
  real ( kind = rk ), parameter :: e5 = 0.0102930D+00
  real ( kind = rk ), parameter :: euler = 2.71828182845904D+00
  real ( kind = rk ) p
  real ( kind = rk ) q
  real ( kind = rk ) q0
  real ( kind = rk ), parameter :: q1 =   0.04166669D+00
  real ( kind = rk ), parameter :: q2 =   0.02083148D+00
  real ( kind = rk ), parameter :: q3 =   0.00801191D+00
  real ( kind = rk ), parameter :: q4 =   0.00144121D+00
  real ( kind = rk ), parameter :: q5 = - 0.00007388D+00
  real ( kind = rk ), parameter :: q6 =   0.00024511D+00
  real ( kind = rk ), parameter :: q7 =   0.00024240D+00
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) si
  real ( kind = rk ) s2
  real ( kind = rk ) t
  real ( kind = rk ) u
  real ( kind = rk ) v
  real ( kind = rk ) w
  real ( kind = rk ) x
!
!  Allow C = 0.
!
  if ( c == 0.0D+00 ) then
    x = a
    return
  end if
!
!  C < 1.
!
  if ( c < 1.0D+00 ) then

    do

      call random_number ( harvest = u )
      t = 1.0D+00 + c / euler
      p = u * t

      call exponential_01_sample ( s )

      if ( p < 1.0D+00 ) then
        x = exp ( log ( p ) / c )
        if ( x <= s ) then
          exit
        end if
      else
        x = - log ( ( t - p ) / c )
        if ( ( 1.0D+00 - c ) * log ( x ) <= s ) then
          exit
        end if
      end if

    end do

    x = a + b * x
    return
!
!  1 <= C.
!
  else

    s2 = c - 0.5D+00
    s = sqrt ( c - 0.5D+00 )
    d = sqrt ( 32.0D+00 ) - 12.0D+00 * sqrt ( c - 0.5D+00 )

    call normal_01_sample ( t )
    x = ( sqrt ( c - 0.5D+00 ) + 0.5D+00 * t )**2

    if ( 0.0D+00 <= t ) then
      x = a + b * x
      return
    end if

    call random_number ( harvest = u )

    if ( d * u <= t**3 ) then
      x = a + b * x
      return
    end if

    r = 1.0D+00 / c

    q0 = ( ( ( ( ( ( &
           q7   * r &
         + q6 ) * r &
         + q5 ) * r &
         + q4 ) * r &
         + q3 ) * r &
         + q2 ) * r &
         + q1 ) * r

    if ( c <= 3.686D+00 ) then
      bcoef = 0.463D+00 + s - 0.178D+00 * s2
      si = 1.235D+00
      co = 0.195D+00 / s - 0.079D+00 + 0.016D+00 * s
    else if ( c <= 13.022D+00 ) then
      bcoef = 1.654D+00 + 0.0076D+00 * s2
      si = 1.68D+00 / s + 0.275D+00
      co = 0.062D+00 / s + 0.024D+00
    else
      bcoef = 1.77D+00
      si = 0.75D+00
      co = 0.1515D+00 / s
    end if

    if ( 0.0D+00 < sqrt ( c - 0.5D+00 ) + 0.5D+00 * t ) then

      v = 0.5D+00 * t / s

      if ( 0.25D+00 < abs ( v ) ) then
        q = q0 - s * t + 0.25D+00 * t * t + 2.0D+00 * s2 * log ( 1.0D+00 + v )
      else
        q = q0 + 0.5D+00 * t**2 * ( ( ( ( ( ( &
               a7   * v &
             + a6 ) * v &
             + a5 ) * v &
             + a4 ) * v &
             + a3 ) * v &
             + a2 ) * v &
             + a1 ) * v
      end if

      if ( log ( 1.0D+00 - u ) <= q ) then
        x = a + b * x
        return
      end if

    end if

    do

      call exponential_01_sample ( e )

      call random_number ( harvest = u )
 
      u = 2.0D+00 * u - 1.0D+00
      t = bcoef + sign ( si * e, u )

      if ( -0.7187449D+00 <= t ) then

        v = 0.5D+00 * t / s

        if ( 0.25D+00 < abs ( v ) ) then
          q = q0 - s * t + 0.25D+00 * t**2 + 2.0D+00 * s2 * log ( 1.0D+00 + v )
        else
          q = q0 + 0.5D+00 * t**2 * ( ( ( ( ( ( &
               a7   * v &
             + a6 ) * v &
             + a5 ) * v &
             + a4 ) * v &
             + a3 ) * v &
             + a2 ) * v &
             + a1 ) * v
        end if

        if ( 0.0D+00 < q ) then

          if ( 0.5D+00 < q ) then
            w = exp ( q ) - 1.0D+00
          else
            w = ( ( ( ( &
                    e5   * q &
                  + e4 ) * q &
                  + e3 ) * q &
                  + e2 ) * q &
                  + e1 ) * q
          end if

          if ( co * abs ( u ) <= w * exp ( e - 0.5D+00 * t**2 ) ) then
            x = a + b * ( s + 0.5D+00 * t )**2
            return
          end if

        end if

      end if

    end do

  end if

  return
end
subroutine i4_next ( line, ival, done )

!*****************************************************************************80
!
!! i4_next() "reads" integers from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) LINE, a string, presumably containing
!    integers.  These may be separated by spaces or commas.
!
!    logical DONE:, on input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!
!  Output:
!
!    integer IVAL.  If DONE is FALSE, then IVAL contains
!    the "next" integer read from LINE.  If DONE is TRUE, then
!    IVAL is zero.
!
!    logical DONE: the routine sets DONE to FALSE if another integer
!    was read, or TRUE if no more integers could be read.
!
  implicit none

  logical done
  integer ierror
  integer ival
  integer lchar
  character ( len = * ) line
  integer, save :: next = 1

  ival = 0

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( len ( line ) < next ) then
    done = .true.
    return
  end if

  call s_to_i4 ( line(next:), ival, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine mixture_read ( acid_num, acid_sym, beta, beta_sum, comp_label, &
  comp_max, comp_num, comp_weight, ierror, iunit )

!*****************************************************************************80
!
!! mixture_read() reads the Dirichlet mixture parameters from a file.
!
!  Discussion:
!
!    The data in the file is delimited by keywords.
!
!    The first lines (not necessarily in order!) may include
!
!      ClassName = string
!      NumDistr = N           the number of components in the mixture.
!      Alphabet = string
!      Order = A C D E ...    the order of the amino acids.
!      AlphaChar = 20
!      NumDistr = 9           the number of distributions
!      EndClassName = string
!
!    For each component, there are four lines:
!
!      Number= N              the component number, starting with 0
!      Mixture= N             the mixture weight, out of a total of 1.0
!      Alpha=  |A| A1 A2 ...  the parameter sum, and individual parameters
!      Comment=               a comment, which describes the frequencies.
!
!    In the comment, the symbol "><" indicates the mean background frequency;
!    residues to the left of that symbol occur more frequently
!    than background, residues to the right less frequently.  Commas separate
!    residues differing in frequency by a factor of 2.
!
!    For example, the comment
!      S A T , C G P >< N V M , Q H R I K F L D W , E Y
!    indicates that for this component, the frequency of
!    proline is just above the mean, and serine, alanine and
!    threonine are twice as frequent in this component than they
!    are on average.  By contrast, tyrosine and glutamic acid are
!    between 4 and 8 times less likely in this component than on
!    average.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer ACID_NUM, the number of amino acids.
!
!    integer COMP_MAX, the maximum number of Dirichlet 
!    mixture components.
!
!    integer IUNIT, the FORTRAN unit from which the 
!    data is to be read.
!
!  Output:
!
!    character ACID_SYM(ACID_NUM), the one letter amino acid codes.
!
!    real ( kind = rk ) BETA(ACID_NUM,COMP_MAX); BETA(I,J) is the 
!    parameter for the J-th acid in the I-th Dirichlet mixture component.
!
!    real ( kind = rk ) BETA_SUM(COMP_MAX), the sum of the values of
!    BETA(ACID_I,COMP_I) for a given component COMP_I.
!
!    integer COMP_LABEL(COMP_NUM), the label of each 
!    component.  Normally, component I has label I.
!
!    integer COMP_NUM, the number of components in the 
!    Dirichlet mixture.
!
!    real ( kind = rk ) COMP_WEIGHT(COMP_NUM), the mixture weight of 
!    each component.  These values should be nonnegative, and sum to 1.  
!    They represent the relative proportion of each component in the mixture.
!
!    integer IERROR, error indicator.
!    0: no error occurred; nonzero: an error occurred.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer acid_num
  integer comp_max

  integer acid_i
  character acid_sym(acid_num)
  real ( kind = rk ) beta(acid_num,comp_max)
  real ( kind = rk ) beta_sum(comp_max)
  integer comp_i
  integer comp_label(comp_max)
  integer comp_num
  real ( kind = rk ) comp_weight(comp_max)
  logical done
  integer iequal
  integer ierror
  integer iunit
  integer ngoofy
  integer nrec
  logical s_begin
  character ( len = 500 ) string

  ierror = 0
  comp_i = 0
  comp_num = 0
  nrec = 0
  ngoofy = 0

10    continue

  read ( iunit, '(a)', end = 20 ) string
  nrec = nrec + 1
!
!  Ignore blank lines.
!
  if ( string == ' ' ) then
!
!  Ignore the CLASSNAME field.
!
  else if ( s_begin ( string, 'CLASSNAME' ) ) then
!
!  Ignore the ENDCLASSNAME field.
!
  else if ( s_begin ( string, 'ENDCLASSNAME' ) ) then
!
!  Ignore the NAME field.
!
  else if ( s_begin ( string, 'NAME' ) ) then
!
!  Ignore the ALPHABET field.
!
  else if ( s_begin ( string, 'ALPHABET' ) ) then
!
!  Read the ORDER field, since it tells us how to interpret the ALPHA's.
!
  else if ( s_begin ( string, 'ORDER' ) ) then

    iequal = index ( string, '=' )
    done = .true.
    do acid_i = 1, acid_num
      call ch_next ( string(iequal+1:), acid_sym(acid_i), done )
    end do
!
!  Ignore the ALPHACHAR field.
!
  else if ( s_begin ( string, 'ALPHACHAR' ) ) then
!
!  Read the NUMDISTR field.
!
  else if ( s_begin ( string, 'NUMDISTR' ) ) then

    iequal = index ( string, '=' )
    done = .true.
    call i4_next ( string(iequal+1:), comp_num, done )

    if ( comp_num < 1 ) then
      ierror = 1
      return
    else if ( comp_max < comp_num ) then
      ierror = 2
      return
    end if
!
!  Read the NUMBER field.
!
  else if ( s_begin ( string, 'NUMBER' ) ) then

    comp_i = comp_i + 1

    if ( comp_num < comp_i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MIXTURE_READ - Fatal error!'
      write ( *, '(a,i6)' ) '  Number of components = ', comp_i
      write ( *, '(a,i6)' ) '  exceeding reported value of ', comp_num
      stop
    end if

    iequal = index ( string, '=' )
    done = .true.
    call i4_next ( string(iequal+1:), comp_label(comp_i), done )
!
!  Read the MIXTURE field.
!
  else if ( s_begin ( string, 'MIXTURE' ) ) then

    iequal = index ( string, '=' )
    done = .true.
    call r8_next ( string(iequal+1:), comp_weight(comp_i), done )
!
!  Read the ALPHA field.
!
  else if ( s_begin ( string, 'ALPHA' ) ) then

    iequal = index ( string, '=' )
    done = .true.
    call r8_next ( string(iequal+1:), beta_sum(comp_i), done )
    do acid_i = 1, acid_num
      call r8_next ( string(iequal+1:), beta(acid_i,comp_i), done )
    end do
!
!  Ignore the COMMENT field.
!
  else if ( s_begin ( string, 'COMMENT' ) ) then
!
!  Unexpected field:
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MIXTURE_READ - Warning!'
    write ( *, '(a)' ) '  Goofy record: '
    write ( *, '(a)' ) trim ( string )

    ngoofy = ngoofy + 1

  end if

  go to 10

20    continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MIXTURE_READ - Note:'
  write ( *, '(a,i6)' ) '  Number of records read was ', nrec
  write ( *, '(a,i6)' ) '  Number of goofy records was ', ngoofy

  return
end
subroutine multinomial_sample ( a, b, c, x )

!*****************************************************************************80
!
!! multinomial_sample() samples the Multinomial PDF.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Luc Devroye,
!    Non-Uniform Random Variate Generation,
!    Springer-Verlag, New York, 1986, page 559.
!
!  Input:
!
!    integer A, the total number of trials.
!    0 <= A.
!
!    integer B, the number of outcomes possible on
!    one trial.  1 <= B.
!
!    real ( kind = rk ) C(B).  C(I) is the probability of outcome I on
!    any trial.
!    0.0 <= C(I) <= 1.0,
!    sum ( 1 <= I <= B) C(I) = 1.0.
!
!  Output:
!
!    integer X(B); X(I) is the number of
!    occurrences of event I during the N trials.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer b

  integer a
  real ( kind = rk ) c(b)
  integer ifactor
  integer ntot
  real ( kind = rk ) prob
  real ( kind = rk ) sum2
  integer x(b)

  ntot = a

  sum2 = 1.0D+00

  x(1:b) = 0

  do ifactor = 1, b - 1

    prob = c(ifactor) / sum2
!
!  Generate a binomial random deviate for NTOT trials with 
!  single trial success probability PROB.
!
    call binomial_sample ( ntot, prob, x(ifactor) )

    ntot = ntot - x(ifactor)
    if ( ntot <= 0 ) then
      return
    end if

    sum2 = sum2 - c(ifactor)

  end do
!
!  The last factor gets what's left.
!
  x(b) = ntot

  return
end
subroutine normal_01_sample ( x )

!*****************************************************************************80
!
!! normal_01_sample() samples the standard normal probability distribution.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has 
!    mean 0 and standard deviation 1.
!
!    The Box-Muller method is used, which is efficient, but 
!    generates two values at a time.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) X, a sample of the standard normal PDF.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  integer, save :: used = -1
  real ( kind = rk ) x
  real ( kind = rk ), save :: y = 0.0D+00

  if ( used == -1 ) then
    used = 0
  end if
!
!  If we've used an even number of values so far, generate two more,
!  return one and save one.
!
  if ( mod ( used, 2 ) == 0 ) then

    call random_number ( harvest = r1 )
    call random_number ( harvest = r2 )
 
    x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  Otherwise, return the second, saved, value.
!
  else

    x = y

  end if

  used = used + 1

  return
end
subroutine r8_next ( line, rval, done )

!*****************************************************************************80
!
!! r8_next() "reads" real numbers from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) LINE, a string, presumably containing real
!    numbers.  These may be separated by spaces or commas.
!
!    logical DONE: on input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!
!  Output:
!
!    real ( kind = rk ) RVAL.  If DONE is FALSE, then RVAL contains the
!    "next" real value read from LINE.  If DONE is TRUE, then
!    RVAL is zero.
!
!    logical DONE: the routine sets DONE to FALSE if another real
!    value was read, or TRUE if no more reals could be read.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical done
  integer ierror
  integer lchar
  character ( len = * ) line
  integer, save :: next = 1
  real ( kind = rk ) rval

  rval = 0.0D+00

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( len ( line ) < next ) then
    done = .true.
    return
  end if

  call s_to_r8 ( line(next:), rval, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! r8vec_print() prints a real vector, with an optional title.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    real ( kind = rk ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r8vec_unit_sum ( n, a )

!*****************************************************************************80
!
!! r8vec_unit_sum() normalizes a real vector to have unit sum.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk ) A(N), the vector to be normalized.  
!
!  Output:
!
!    real ( kind = rk ) A(N), the entries of A should have unit sum.  
!    However, if the input vector has zero sum, the routine halts.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_sum

  a_sum = sum ( a(1:n) )

  if ( a_sum == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIT_SUM - Fatal error!'
    write ( *, '(a)' ) '  The vector entries sum to 0.'
    stop
  end if

  a(1:n) = a(1:n) / a_sum

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! r8_gamma_log() calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Discussion:
!
!    The program uses rational functions that theoretically approximate
!    log ( GAMMA(X) ) to at least 18 significant decimal digits.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Mathematics of Computation,
!    Volume 21, 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    Hart, Ward Cheney, Charles Lawson, Maehly, Charles Mesztenyi, 
!    John Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Input:
!
!    real ( kind = rk ) X, the argument of the Gamma function. 
!    X must be positive.
!
!  Output:
!
!    real ( kind = rk ) R8_GAMMA_LOG, the logarithm of the Gamma 
!    function of X.  If X <= 0.0, or if overflow would occur, the
!    program returns the value HUGE().
!
!  Machine-dependent constants:
!
!    BETA   - radix for the floating-point representation.
!
!    MAXEXP - the smallest positive power of BETA that overflows.
!
!    XBIG   - largest argument for which LN(GAMMA(X)) is representable
!             in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA^MAXEXP.
!
!    XINF   - largest machine representable floating-point number;
!             approximately BETA^MAXEXP.
!
!    FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!    Approximate values for some important machines are:
!
!                              BETA      MAXEXP         XBIG
!
!    CRAY-1        (S.P.)        2        8191       9.62D+2461
!    Cyber 180/855
!      under NOS   (S.P.)        2        1070       1.72D+319
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)        2         128       4.08D+36
!    IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!    IBM 3033      (D.P.)       16          63       4.29D+73
!    VAX D-Format  (D.P.)        2         127       2.05D+36
!    VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            FRTBIG
!
!    CRAY-1        (S.P.)   3.13D+615
!    Cyber 180/855
!      under NOS   (S.P.)   6.44D+79
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)   1.42D+9
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)   2.25D+76
!    IBM 3033      (D.P.)   2.56D+18
!    VAX D-Format  (D.P.)   1.20D+9
!    VAX G-Format  (D.P.)   1.89D+76
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = rk ) corr
  real ( kind = rk ), parameter :: d1 = - 5.772156649015328605195174D-01
  real ( kind = rk ), parameter :: d2 =   4.227843350984671393993777D-01
  real ( kind = rk ), parameter :: d4 =   1.791759469228055000094023D+00
  real ( kind = rk ) eps
  real ( kind = rk ), parameter :: frtbig = 1.42D+09
  integer i
  real ( kind = rk ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = rk ), parameter :: pnt68 = 0.6796875D+00
  real ( kind = rk ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = rk ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = rk ) r8_gamma_log
  real ( kind = rk ) res
  real ( kind = rk ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = rk ) x
  real ( kind = rk ), parameter :: xbig = 4.08D+36
  real ( kind = rk ) xden
  real ( kind = rk ) xm1
  real ( kind = rk ) xm2
  real ( kind = rk ) xm4
  real ( kind = rk ) xnum
  real ( kind = rk ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    r8_gamma_log = huge ( r8_gamma_log )
    return
  end if

  eps = epsilon ( eps )

  if ( x <= eps ) then

    res = - log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0D+00
      xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if

    if ( x <= 0.5D+00 .or. pnt68 <= x ) then

      xden = 1.0D+00
      xnum = 0.0D+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0D+00 ) then

    xm2 = x - 2.0D+00
    xden = 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0D+00 ) then

    xm4 = x - 4.0D+00
    xden = -1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0D+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )

  end if

  r8_gamma_log = res

  return
end
function s_begin ( s1, s2 )

!*****************************************************************************80
!
!! s_begin() is TRUE if one string matches the beginning of the other.
!
!  Discussion:
!
!    The strings are compared, ignoring blanks, spaces and capitalization.
!
!  Example:
!
!     S1              S2      S_BEGIN
!
!    'Bob'          'BOB'     TRUE
!    '  B  o b '    ' bo b'   TRUE
!    'Bob'          'Bobby'   TRUE
!    'Bobo'         'Bobb'    FALSE
!    ' '            'Bob'     FALSE    (Do not allow a blank to match
!                                       anything but another blank string.)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) S1, S2, the strings to be compared.
!
!  Output:
!
!    logical S_BEGIN, is TRUE if the strings match up to
!    the end of the shorter string, ignoring case.
!
  implicit none

  logical ch_eqi
  integer i1
  integer i2
  logical s_begin
  character ( len = * ) s1
  integer s1_length
  character ( len = * ) s2
  integer s2_length

  s1_length = len_trim ( s1 )
  s2_length = len_trim ( s2 )
!
!  If either string is blank, then both must be blank to match.
!  Otherwise, a blank string matches anything, which is not
!  what most people want.
!
  if ( s1_length == 0 .or. s2_length == 0 ) then

    if ( s1_length == 0 .and. s2_length == 0 ) then
      s_begin = .true.
    else
      s_begin = .false.
    end if

    return

  end if

  i1 = 0
  i2 = 0
!
!  Find the next nonblank in S1.
!
  do

    do

      i1 = i1 + 1

      if ( s1_length < i1 ) then
        s_begin = .true.
        return
      end if

      if ( s1(i1:i1) /= ' ' ) then
        exit
      end if

    end do
!
!  Find the next nonblank in S2.
!
    do

      i2 = i2 + 1

      if ( s2_length < i2 ) then
        s_begin = .true.
        return
      end if

      if ( s2(i2:i2) /= ' ' ) then
        exit
      end if

    end do
!
!  If the characters match, get the next pair.
!
    if ( .not. ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
      exit
    end if

  end do

  s_begin = .false.

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! s_to_i4() reads an I4 from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) S, a string to be examined.
!
!  Output:
!
!    integer VALUE, the value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    integer LENGTH, the number of characters
!    of S used to make the integer.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer length
  character ( len = * ) s
  integer state
  character :: TAB = achar ( 9 )
  integer value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! s_to_r8() reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = rk )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!  Output:
!
!    real ( kind = rk ) DVAL, the value read from the string.
!
!    integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    integer LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character c
  logical ch_eqi
  real ( kind = rk ) dval
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer ndig
  real ( kind = rk ) rbot
  real ( kind = rk ) rexp
  real ( kind = rk ) rtop
  character ( len = * ) s
  integer s_length
  character :: TAB = achar ( 9 )

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = rk )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = rk ) &
        / real ( jbot, kind = rk ) )
    end if
  end if

  dval = real ( isgn, kind = rk ) * rexp * rtop / rbot

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
 
