cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Followinf subroutines are copied from Binbin Zhou Grid2grid program
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	SUBROUTINE ST_RMBL  ( string, outstr, length, iret )
C************************************************************************
C* ST_RMBL								*
C*									*
C* This subroutine removes spaces and tabs from a string.  The input	*
C* and output strings may be the same variable.				*
C*									*
C* ST_RMBL  ( STRING, OUTSTR, LENGTH, IRET )				*
C*									*
C* Input parameters:							*
C*	STRING		CHAR*		String				*
C*									*
C* Output parameters:							*
C*	OUTSTR		CHAR*		String without blanks		*
C*	LENGTH		INTEGER		Length of output string		*
C*	IRET		INTEGER		Return code			*
C*				   	 0 = normal return 		*
C**									*
C* Log:									*
C* M. desJardins/GSFC	 8/84						*
C* M. desJardins/GSFC	 6/88	Documentation				*
C* M. desJardins/NMC	 3/92	Add temporary variable ttt		*
C* L. Sager/NCEP         2/96   Increased size of sss and ttt           *
C* D. Kidwell/NCEP      10/96   Ported to Cray                          *
C* B. Zhou              05/2005 Adapted to Grid2Grid                    * 
C************************************************************************
C* GEMPRM.PRM
C*
C* This include file contains parameter definitions for the GEMPAK
C* software routines in the ST_ and PR_ libraries.
C*
C* CRAY version
C**
C* Log:  
C*	Kidwell/NCEP	07/96	Adapted a subset of gemprm.prm for Cray
C************************************************************************
C!
C!	Missing data definitions
C!
 	PARAMETER	( RMISSD = -9999.0 )
C!						Missing data value
	PARAMETER	( RDIFFD =  0.1    )
C!						Missing value fuzziness
	PARAMETER	( IMISSD = -9999   )
C!						Missing integer value
	LOGICAL		  ERMISS
C!						Declare for stmt func
C!
C! 	Physical and mathematical constants
C!
	PARAMETER       ( PI = 3.14159265  )
C!                                              PI
	PARAMETER       ( DTR = PI / 180.  )
	PARAMETER       ( RTD = 180. / PI  )
C!                                              Degrees <--> Radians
	PARAMETER	( GRAVTY = 9.80616  )
C!						Acceleration of gravity
	PARAMETER	( RDGAS  = 287.04   )
	PARAMETER	( RKAP   = RDGAS / GRAVTY )
C!						Gas constant of dry air
	PARAMETER	( RKAPPA = 2. / 7. )
C!						Poisson constant
	PARAMETER	( GAMUSD = 6.5 )
C!						US std atmos lapse rate
	PARAMETER	( TMCK   = 273.15 )
C!						Centigrade -> Kelvin
C!
C!	ASCII character constants 
C!
C!	Since the Cray does not allow the use of a function (e.g.,
C!	CHAR) to define a parameter, nor does it allow a character
C!	to be defined directly as a hex (X) value, the convolutions
C!      below are necessary to define the character values for Cray.
C!
C
	CHARACTER * 1 chnull, chtab, chspac, chtlda
C
	CHARACTER * 8 c8null, c8tab, c8spac, c8tlda
C
	INTEGER       iigemc ( 4 )
C
	EQUIVALENCE ( chnull, c8null (8:8) ), ( chtab , c8tab  (8:8) ),
     +              ( chspac, c8spac (8:8) ), ( chtlda, c8tlda (8:8) ) 
C
	EQUIVALENCE ( iigemc ( 1), c8null ), ( iigemc ( 2), c8tab  ),
     +              ( iigemc ( 3), c8spac ), ( iigemc ( 4), c8tlda )
C
	DATA iigemc / X'00',    X'09',    X'20',    X'7E' /
C                     Null      Tab       Space     Tilda
C!
C*
	CHARACTER*(*)	string, outstr
C*
	CHARACTER	c*1, sss*160, ttt*160
C-----------------------------------------------------------------------
	iret   = 0
	length = 0
	sss    = string
	ttt    = ' '
C
C*	Get length of input string.
C
	CALL ST_LSTR  ( sss, lens, iret )
C
C*	Check each character to see if it is a blank.
C
	DO  i = 1, lens
	    c = sss (i:i)
	    IF  ( ( c .ne. CHSPAC ) .and. ( c .ne. CHTAB ) )  THEN
		length = length + 1
		ttt ( length : length ) = c
	    END IF
	END DO
C*
	outstr = ttt
C*
	RETURN
	END


	SUBROUTINE ST_NUMB  ( string, ival, iret )
C************************************************************************
C* ST_NUMB                                       			*
C*									*
C* This subroutine converts a string into an integer.			* 
C*									*
C* ST_NUMB  ( STRING, IVAL, IRET )					*
C*									*
C* Input parameters:							*
C*	STRING		CHAR*		String				*
C*									*
C* Output parameters:							*
C*	IVAL		INTEGER		Integer value      		*
C*	IRET		INTEGER		Return code			*
C*				 	  0 = normal return		*
C*					 -2 = conversion error		*
C** 									*
C* Log:									*
C* M. desJardins/NMC	 3/92	Rewritten to avoid special cases	*
C* D. Kidwell/NCEP      10/96   Ported to Cray                          *
C* B. Zhou              05/2005 Adapted to Grid2Grid                    *
C************************************************************************
C* GEMPRM.PRM
C*
C* This include file contains parameter definitions for the GEMPAK
C* software routines in the ST_ and PR_ libraries.
C*
C* CRAY version
C**
C* Log:  
C*	Kidwell/NCEP	07/96	Adapted a subset of gemprm.prm for Cray
C************************************************************************
C!
C!	Missing data definitions
C!
 	PARAMETER	( RMISSD = -9999.0 )
C!						Missing data value
	PARAMETER	( RDIFFD =  0.1    )
C!						Missing value fuzziness
	PARAMETER	( IMISSD = -9999   )
C!						Missing integer value
	LOGICAL		  ERMISS
C!						Declare for stmt func
C!
C! 	Physical and mathematical constants
C!
	PARAMETER       ( PI = 3.14159265  )
C!                                              PI
	PARAMETER       ( DTR = PI / 180.  )
	PARAMETER       ( RTD = 180. / PI  )
C!                                              Degrees <--> Radians
	PARAMETER	( GRAVTY = 9.80616  )
C!						Acceleration of gravity
	PARAMETER	( RDGAS  = 287.04   )
	PARAMETER	( RKAP   = RDGAS / GRAVTY )
C!						Gas constant of dry air
	PARAMETER	( RKAPPA = 2. / 7. )
C!						Poisson constant
	PARAMETER	( GAMUSD = 6.5 )
C!						US std atmos lapse rate
	PARAMETER	( TMCK   = 273.15 )
C!						Centigrade -> Kelvin
C!
C!	ASCII character constants 
C!
C!	Since the Cray does not allow the use of a function (e.g.,
C!	CHAR) to define a parameter, nor does it allow a character
C!	to be defined directly as a hex (X) value, the convolutions
C!      below are necessary to define the character values for Cray.
C!
C
	CHARACTER * 1 chnull, chtab, chspac, chtlda
C
	CHARACTER * 8 c8null, c8tab, c8spac, c8tlda
C
	INTEGER       iigemc ( 4 )
C
	EQUIVALENCE ( chnull, c8null (8:8) ), ( chtab , c8tab  (8:8) ),
     +              ( chspac, c8spac (8:8) ), ( chtlda, c8tlda (8:8) ) 
C
	EQUIVALENCE ( iigemc ( 1), c8null ), ( iigemc ( 2), c8tab  ),
     +              ( iigemc ( 3), c8spac ), ( iigemc ( 4), c8tlda )
C
	DATA iigemc / X'00',    X'09',    X'20',    X'7E' /
C                     Null      Tab       Space     Tilda
C!
C*
	CHARACTER*(*) 	string
C*
	CHARACTER	sss*12, c*1
	LOGICAL		good, plus
C------------------------------------------------------------------------
	iret = -2
	ival = IMISSD
C
C*	Remove blanks from string.
C
CZeus	CALL ST_RMBL  ( string, sss, lens, ier )
        lens=len_trim(string)
        sss=string(1:lens)

C
C*	Check for empty string.
C
	IF  ( lens .eq. 0 )  RETURN
C
C*	If last character is period, remove it.
C
	IF  ( sss (lens:lens) .eq. '.' )  THEN
	    sss (lens:lens) = ' '
	    lens = lens - 1
	END IF
	IF  ( lens .eq. 0 )  RETURN
C
C*	Check for + or - in first character.
C
	IF  ( sss (1:1) .eq. '+' )  THEN
	    ibeg = 2
	    plus = .true.
	  ELSE IF  ( sss (1:1) .eq. '-' )  THEN
	    ibeg = 2
	    plus = .false.
	  ELSE
	    ibeg = 1
	    plus = .true.
	END IF
	IF  ( ibeg .gt. lens )  RETURN
C
C*	Now loop through all characters and turn into integer.
C
	ival0 = ICHAR ( '0' )
	ival  = 0
	good  = .true.
	i     = ibeg
	DO  WHILE  ( good .and. ( i .le. lens ) )
	    c     = sss (i:i)
	    ivalc = ICHAR ( c ) - ival0
	    IF  ( ( ivalc .ge. 0 ) .and. ( ivalc .le. 9 ) )  THEN
		ival  = ival * 10 + ivalc
		i     = i + 1
	      ELSE
		good  = .false.
	    END IF
	END DO
C
C*	Check for good value and add sign.
C
	IF  ( .not. good )  THEN
	    ival = IMISSD
	  ELSE IF  ( plus )  THEN
	    iret = 0
	  ELSE
	    iret = 0
	    ival = -ival
	END IF
C*
	RETURN
	END


	SUBROUTINE ST_CLST  ( string, sep, cdef, nexp, carr, num, iret )
C************************************************************************
C* ST_CLST								*
C*									*
C* This subroutine breaks a string containing a list of strings into	*
C* an array of strings.  The separator for the strings is input as SEP.	*
C* If the separator is a blank, multiple blanks will be changed to	*
C* single blanks before the string is processed.  If null strings	*
C* are encountered or fewer than NEXP strings are found in the		*
C* string, the appropriate CARR locations are set to CDEF.		*
C*									*
C* ST_CLST  ( STRING, SEP, CDEF, NEXP, CARR, NUM, IRET )		*
C*									*
C* Input parameters:							*
C*	STRING		CHAR*		String 				*
C*	SEP		CHAR*1		Separator			*
C*	CDEF		CHAR*		Default string 			*
C*	NEXP		INTEGER 	Number of expected values 	*
C*									*
C* Output parameters:							*
C*	CARR  (NUM)	CHAR*		Array of strings 		*
C*	NUM		INTEGER 	Number of strings returned	*
C*	IRET		INTEGER 	Return code			*
C*				   	  1 = more than NEXP values	*
C*				   	  0 = normal return		*
C**									*
C* Log:									*
C* M. Goodman/RDS	10/84	Original source for STLIST		*
C* M. desJardins/GSFC	11/84						*
C* M. desJardins/GSFC	 2/85	Modified for ST_CLST			*
C* M. desJardins/GSFC	 6/88	Documentation				*
C* L. Sager/NCEP         2/96   Increased size of strbuf                *
C* D. Kidwell/NCEP      10/96   Ported to Cray                          *
C* B. Zhou              05/2005 Adapted to Grid2Grid                    *
C************************************************************************
	CHARACTER*(*) 	string, sep, carr (*), cdef
C*
	CHARACTER	strbuf*160, cchar*1
C------------------------------------------------------------------------
	iret = 0
	num  = 0
C
C*	Remove blanks from the input string if the separator is not
C*	a blank.
C
	IF  ( sep .ne. ' ' )  THEN
	    CALL ST_RMBL  ( string, strbuf, isize, iret )
	  ELSE
	    CALL ST_RXBL  ( string, strbuf, isize, iret )
	END IF
C
C*	Initialize output array.
C
	DO  i = 1, nexp
	    carr (i) = cdef
	END DO
C
C*	Check for zero length input string.
C
	IF  ( isize .eq. 0 )  THEN
	    num = 0
C
C*	    Check for separator and find list elements.
C
	  ELSE
	    cchar  = sep
	    iend   = 0
	    ibegin = 1
	    DO WHILE  ( ibegin .le. isize )
	        loc = INDEX  ( strbuf ( ibegin: ), cchar )
	        IF  ( loc .eq. 0 )  THEN
		    iend = isize + 1
	          ELSE
		    iend = ibegin + loc - 1
	        END IF
C
C*	        Add into output list.  Check that num <= nexp.
C
	        IF  ( num .ge. nexp )  THEN
		    iret = 1
		  ELSE
		    num = num + 1
	    	    IF  ( ibegin .ne. iend )  THEN
			carr ( num ) = strbuf ( ibegin : iend-1 )
	    	    END IF
		END IF

CZeus modify to skip space: by Binbin Zhoucccccccccccccccccccccccccc
! mptest to prevent bounds exceedance
	if (iend+1 .le. 160)  then
         do while (strbuf(iend+1:iend+1) .eq. ' ' )
           iend = iend + 1
           if (iend .eq. 160) exit
         end do
	end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccc


	    ibegin = iend + 1
	    END DO
	END IF
C*
	RETURN
	END


	SUBROUTINE ST_LDSP  ( string, outstr, ncout, iret )
C************************************************************************
C* ST_LDSP								*
C*									*
C* This subroutine deletes the leading spaces and tabs in a string.	*
C* The input and output strings may be the same variable.		*
C*									*
C* ST_LDSP  ( STRING, OUTSTR, NCOUT, IRET )				*
C*									*
C* Input parameters:							*
C*	STRING		CHAR*		String				*
C*									*
C* Output parameters:							*
C*	OUTSTR		CHAR*		Output string			*
C*	NCOUT		INTEGER		Number of characters output	*
C*	IRET		INTEGER		Return code			*
C*					  0 = normal return		*
C**									*
C* Log:									*
C* I. Graffman/RDS	 2/84	Use new GEMPAK routines			*
C* M. desJardins/GSFC	11/84	Fixed					*
C* M. desJardins/GSFC	 6/88	Documentation				*
C* L. Sager/NCEP         2/96   Increased size of stbuf                 *
C* D. Kidwell/NCEP      10/96   Ported to Cray                          *
C* B. Zhou              05/2005 Adapted to Grid2Grid                    *
C************************************************************************
C* GEMPRM.PRM
C*
C* This include file contains parameter definitions for the GEMPAK
C* software routines in the ST_ and PR_ libraries.
C*
C* CRAY version
C**
C* Log:  
C*	Kidwell/NCEP	07/96	Adapted a subset of gemprm.prm for Cray
C************************************************************************
C!
C!	Missing data definitions
C!
 	PARAMETER	( RMISSD = -9999.0 )
C!						Missing data value
	PARAMETER	( RDIFFD =  0.1    )
C!						Missing value fuzziness
	PARAMETER	( IMISSD = -9999   )
C!						Missing integer value
	LOGICAL		  ERMISS
C!						Declare for stmt func
C!
C! 	Physical and mathematical constants
C!
	PARAMETER       ( PI = 3.14159265  )
C!                                              PI
	PARAMETER       ( DTR = PI / 180.  )
	PARAMETER       ( RTD = 180. / PI  )
C!                                              Degrees <--> Radians
	PARAMETER	( GRAVTY = 9.80616  )
C!						Acceleration of gravity
	PARAMETER	( RDGAS  = 287.04   )
	PARAMETER	( RKAP   = RDGAS / GRAVTY )
C!						Gas constant of dry air
	PARAMETER	( RKAPPA = 2. / 7. )
C!						Poisson constant
	PARAMETER	( GAMUSD = 6.5 )
C!						US std atmos lapse rate
	PARAMETER	( TMCK   = 273.15 )
C!						Centigrade -> Kelvin
C!
C!	ASCII character constants 
C!
C!	Since the Cray does not allow the use of a function (e.g.,
C!	CHAR) to define a parameter, nor does it allow a character
C!	to be defined directly as a hex (X) value, the convolutions
C!      below are necessary to define the character values for Cray.
C!
C
	CHARACTER * 1 chnull, chtab, chspac, chtlda
C
	CHARACTER * 8 c8null, c8tab, c8spac, c8tlda
C
	INTEGER       iigemc ( 4 )
C
	EQUIVALENCE ( chnull, c8null (8:8) ), ( chtab , c8tab  (8:8) ),
     +              ( chspac, c8spac (8:8) ), ( chtlda, c8tlda (8:8) ) 
C
	EQUIVALENCE ( iigemc ( 1), c8null ), ( iigemc ( 2), c8tab  ),
     +              ( iigemc ( 3), c8spac ), ( iigemc ( 4), c8tlda )
C
	DATA iigemc / X'00',    X'09',    X'20',    X'7E' /
C                     Null      Tab       Space     Tilda
C!
C*
	CHARACTER*(*)	string, outstr
C*
	CHARACTER	stbuf*160, c*1
C*-------------------------------------------------------------------------
	stbuf = string
	iret  = 0
C
C*	Get length of string.
C
	CALL ST_LSTR  ( stbuf, lens, iret )
C
C*	If length is non-zero, find first non space.
C
	IF  ( lens .eq. 0 )  THEN
	    ncout  = 0
	    outstr = ' '
	  ELSE
	    jp = 1
	    c  = stbuf ( jp:jp )
C
	    DO WHILE  ( ( ( c .eq. CHSPAC ) .or. ( c .eq. CHTAB ) .or.
     +			  ( c .eq. CHNULL ) ) .and. ( jp .le. lens ) )
		jp = jp + 1
		IF  ( jp .le. lens )  c = stbuf ( jp:jp )
	    ENDDO
C
C*	    Compute length and fill output string.
C
	    ncout = lens - jp + 1
	    IF  ( ncout .gt. 0 )  THEN
		outstr = stbuf ( jp : lens )
	      ELSE
		outstr = ' '
	    END IF
	ENDIF
C*  
	RETURN
	END


	SUBROUTINE ST_LSTR  ( string, lens, iret )
C************************************************************************
C* ST_LSTR								*
C*									*
C* This subroutine returns the number of characters in a string 	*
C* disregarding trailing null characters, tabs and spaces.		*
C*									*
C* ST_LSTR  ( STRING, LENS, IRET )						*
C*									*
C* Input parameters:							*
C*	STRING		CHAR*		String 				*
C*									*
C* Output parameters:							*
C*	LENS		INTEGER 	Length of string		*
C*	IRET		INTEGER		Return code			*
C*				 	 0 = normal return 		*
C**									*
C* Log:									*
C* J. Woytek/GSFC	 6/82 	STR_LNSTR				*
C* I. Graffman/RDS	 2/84 	Fix zero length string handling		*
C* M. desJardins/GSFC	 6/88	Rewritten				*
C* D. Kidwell/NCEP      10/96   Ported to Cray                          *
C* B. Zhou              05/2005 Adapted to Grid2Grid                    *
C************************************************************************
C* GEMPRM.PRM
C*
C* This include file contains parameter definitions for the GEMPAK
C* software routines in the ST_ and PR_ libraries.
C*
C* CRAY version
C**
C* Log:  
C*	Kidwell/NCEP	07/96	Adapted a subset of gemprm.prm for Cray
C************************************************************************
C!
C!	Missing data definitions
C!
 	PARAMETER	( RMISSD = -9999.0 )
C!						Missing data value
	PARAMETER	( RDIFFD =  0.1    )
C!						Missing value fuzziness
	PARAMETER	( IMISSD = -9999   )
C!						Missing integer value
	LOGICAL		  ERMISS
C!						Declare for stmt func
C!
C! 	Physical and mathematical constants
C!
	PARAMETER       ( PI = 3.14159265  )
C!                                              PI
	PARAMETER       ( DTR = PI / 180.  )
	PARAMETER       ( RTD = 180. / PI  )
C!                                              Degrees <--> Radians
	PARAMETER	( GRAVTY = 9.80616  )
C!						Acceleration of gravity
	PARAMETER	( RDGAS  = 287.04   )
	PARAMETER	( RKAP   = RDGAS / GRAVTY )
C!						Gas constant of dry air
	PARAMETER	( RKAPPA = 2. / 7. )
C!						Poisson constant
	PARAMETER	( GAMUSD = 6.5 )
C!						US std atmos lapse rate
	PARAMETER	( TMCK   = 273.15 )
C!						Centigrade -> Kelvin
C!
C!	ASCII character constants 
C!
C!	Since the Cray does not allow the use of a function (e.g.,
C!	CHAR) to define a parameter, nor does it allow a character
C!	to be defined directly as a hex (X) value, the convolutions
C!      below are necessary to define the character values for Cray.
C!
C
	CHARACTER * 1 chnull, chtab, chspac, chtlda
C
	CHARACTER * 8 c8null, c8tab, c8spac, c8tlda
C
	INTEGER       iigemc ( 4 )
C
	EQUIVALENCE ( chnull, c8null (8:8) ), ( chtab , c8tab  (8:8) ),
     +              ( chspac, c8spac (8:8) ), ( chtlda, c8tlda (8:8) ) 
C
	EQUIVALENCE ( iigemc ( 1), c8null ), ( iigemc ( 2), c8tab  ),
     +              ( iigemc ( 3), c8spac ), ( iigemc ( 4), c8tlda )
C
	DATA iigemc / X'00',    X'09',    X'20',    X'7E' /
C                     Null      Tab       Space     Tilda
C!
C*
	CHARACTER*(*)	string
C*
	CHARACTER*1	c
C*------------------------------------------------------------------------
	lens = 0
	iret = 0
C
C*	Get the actual length of the string.
C
	lens = LEN  ( string )
	IF  ( lens .eq. 0 )  RETURN
C
C*	Start at last character and loop backwards.
C
	ip = lens
	DO WHILE  ( ip .gt. 0 )
C
C*	    Get current value of string and check for space, null, tab.
c
	    c = string ( ip : ip )
	    IF  ( ( c .eq. CHSPAC ) .or. ( c .eq. CHNULL ) .or.
     +		  ( c .eq. CHTAB  ) )  THEN
		lens = lens - 1
		ip   = ip - 1
	      ELSE
		ip   = 0
	    END IF
	END DO
C*
	RETURN
	END


	SUBROUTINE ST_RXBL  ( string, outstr, length, iret )
C************************************************************************
C* ST_RXBL								*
C*									*
C* This subroutine removes extra spaces and tabs from a string.  Only	*
C* single blanks will separate substrings.  The input and output 	*
C* strings may be the same variable.					*
C*									*
C* ST_RXBL  ( STRING, OUTSTR, LENGTH, IRET )				*
C*									*
C* Input parameters:							*
C*	STRING		CHAR*		String				*
C*									*
C* Output parameters:							*
C*	OUTSTR		CHAR*		String without blanks		*
C*	LENGTH		INTEGER		Length of output string		*
C*	IRET		INTEGER		Return code			*
C*				   	 0 = normal return 		*
C**									*
C* Log:									*
C* M. desJardins/GSFC	 9/88						*
C* D. Kidwell/NCEP      10/96   Ported to Cray                          *
C* B. Zhou              05/2005 Adapted to Grid2Grid                    *
C************************************************************************
C* GEMPRM.PRM
C*
C* This include file contains parameter definitions for the GEMPAK
C* software routines in the ST_ and PR_ libraries.
C*
C* CRAY version
C**
C* Log:  
C*	Kidwell/NCEP	07/96	Adapted a subset of gemprm.prm for Cray
C************************************************************************
C!
C!	Missing data definitions
C!
 	PARAMETER	( RMISSD = -9999.0 )
C!						Missing data value
	PARAMETER	( RDIFFD =  0.1    )
C!						Missing value fuzziness
	PARAMETER	( IMISSD = -9999   )
C!						Missing integer value
	LOGICAL		  ERMISS
C!						Declare for stmt func
C!
C! 	Physical and mathematical constants
C!
	PARAMETER       ( PI = 3.14159265  )
C!                                              PI
	PARAMETER       ( DTR = PI / 180.  )
	PARAMETER       ( RTD = 180. / PI  )
C!                                              Degrees <--> Radians
	PARAMETER	( GRAVTY = 9.80616  )
C!						Acceleration of gravity
	PARAMETER	( RDGAS  = 287.04   )
	PARAMETER	( RKAP   = RDGAS / GRAVTY )
C!						Gas constant of dry air
	PARAMETER	( RKAPPA = 2. / 7. )
C!						Poisson constant
	PARAMETER	( GAMUSD = 6.5 )
C!						US std atmos lapse rate
	PARAMETER	( TMCK   = 273.15 )
C!						Centigrade -> Kelvin
C!
C!	ASCII character constants 
C!
C!	Since the Cray does not allow the use of a function (e.g.,
C!	CHAR) to define a parameter, nor does it allow a character
C!	to be defined directly as a hex (X) value, the convolutions
C!      below are necessary to define the character values for Cray.
C!
C
	CHARACTER * 1 chnull, chtab, chspac, chtlda
C
	CHARACTER * 8 c8null, c8tab, c8spac, c8tlda
C
	INTEGER       iigemc ( 4 )
C
	EQUIVALENCE ( chnull, c8null (8:8) ), ( chtab , c8tab  (8:8) ),
     +              ( chspac, c8spac (8:8) ), ( chtlda, c8tlda (8:8) ) 
C
	EQUIVALENCE ( iigemc ( 1), c8null ), ( iigemc ( 2), c8tab  ),
     +              ( iigemc ( 3), c8spac ), ( iigemc ( 4), c8tlda )
C
	DATA iigemc / X'00',    X'09',    X'20',    X'7E' /
C                     Null      Tab       Space     Tilda
C!
C*
	CHARACTER*(*)	string, outstr
C-----------------------------------------------------------------------
	length = 0
	iret   = 0
C
C*	Remove leading spaces and tabs.
C
	CALL ST_LDSP  ( string, outstr, isiz, iret )
	IF  ( isiz .le. 0 )  RETURN
C
C*	Remove extra spaces.
C
	ispac  = 0
	length = 0
	DO  j = 1, isiz
	    IF  ( ( outstr (j:j) .ne. CHSPAC )  .and. 
     +		  ( outstr (j:j) .ne. CHTAB  ) )  THEN
	        length = length + 1
	        outstr (length:length) = outstr (j:j)
	        ispac = 0
	      ELSE
	        IF  ( ispac .eq. 0 )  THEN
	            length = length + 1
	            outstr (length:length) = ' '
	            ispac = 1
	        END IF
	    END IF
	END DO
C
C*	Make sure the end of the string is blank.
C
	lens = LEN ( outstr )
	IF  ( lens .gt. length )  outstr ( length+1 : ) = ' '
C*
	RETURN
	END


