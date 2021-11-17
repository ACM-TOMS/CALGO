!
! netCDF VERSION 3 FORTRAN INTERFACE:
!

!
! External netCDF data types:
!
      INTEGER NF_BYTE
      INTEGER NF_INT1
      INTEGER NF_CHAR
      INTEGER NF_SHORT
      INTEGER NF_INT2
      INTEGER NF_INT
      INTEGER NF_FLOAT
      INTEGER NF_REAL
      INTEGER NF_DOUBLE

      PARAMETER (NF_BYTE = 1)
      PARAMETER (NF_INT1 = NF_BYTE)
      PARAMETER (NF_CHAR = 2)
      PARAMETER (NF_SHORT = 3)
      PARAMETER (NF_INT2 = NF_SHORT)
      PARAMETER (NF_INT = 4)
      PARAMETER (NF_FLOAT = 5)
      PARAMETER (NF_REAL = NF_FLOAT)
      PARAMETER (NF_DOUBLE = 6)

!
! Default fill values:
!
      INTEGER		NF_FILL_BYTE
      INTEGER		NF_FILL_INT1
      INTEGER		NF_FILL_CHAR
      INTEGER		NF_FILL_SHORT
      INTEGER		NF_FILL_INT2
      INTEGER		NF_FILL_INT
      REAL		NF_FILL_FLOAT
      REAL		NF_FILL_REAL
      DOUBLEPRECISION	NF_FILL_DOUBLE

      PARAMETER (NF_FILL_BYTE = -127)
      PARAMETER (NF_FILL_INT1 = NF_FILL_BYTE)
      PARAMETER (NF_FILL_CHAR = 0)
      PARAMETER (NF_FILL_SHORT = -32767)
      PARAMETER (NF_FILL_INT2 = NF_FILL_SHORT)
      PARAMETER (NF_FILL_INT = -2147483647)
      PARAMETER (NF_FILL_FLOAT = 9.9692099683868690e+36)
      PARAMETER (NF_FILL_REAL = NF_FILL_FLOAT)
      PARAMETER (NF_FILL_DOUBLE = 9.9692099683868690e+36)

!
! Mode flags for opening and creating a netCDF dataset:
!
      INTEGER NF_NOWRITE
      INTEGER NF_WRITE
      INTEGER NF_CLOBBER
      INTEGER NF_NOCLOBBER
      INTEGER NF_FILL
      INTEGER NF_NOFILL
      INTEGER NF_LOCK
      INTEGER NF_SHARE

      PARAMETER (NF_NOWRITE = 0)
      PARAMETER (NF_WRITE = 1)
      PARAMETER (NF_CLOBBER = 0)
      PARAMETER (NF_NOCLOBBER = 4)
      PARAMETER (NF_FILL = 0)
      PARAMETER (NF_NOFILL = 256)
      PARAMETER (NF_LOCK = 1024)
      PARAMETER (NF_SHARE = 2048)

!
! Size argument for defining an unlimited dimension:
!
      INTEGER NF_UNLIMITED
      PARAMETER (NF_UNLIMITED = 0)

!
! Global Attribute ID:
!
      INTEGER NF_GLOBAL
      PARAMETER (NF_GLOBAL = 0)

!
! Implementation limits:
!
      INTEGER NF_MAX_DIMS
      INTEGER NF_MAX_ATTRS
      INTEGER NF_MAX_VARS
      INTEGER NF_MAX_NAME
      INTEGER NF_MAX_VAR_DIMS

      PARAMETER (NF_MAX_DIMS = 100)
      PARAMETER (NF_MAX_ATTRS = 2000)
      PARAMETER (NF_MAX_VARS = 2000)
      PARAMETER (NF_MAX_NAME = 128)
      PARAMETER (NF_MAX_VAR_DIMS = NF_MAX_DIMS)

!
! Error codes:
!
      INTEGER NF_NOERR
      INTEGER NF_EBADID
      INTEGER NF_EEXIST
      INTEGER NF_EINVAL
      INTEGER NF_EPERM
      INTEGER NF_ENOTINDEFINE
      INTEGER NF_EINDEFINE
      INTEGER NF_EINVALCOORDS
      INTEGER NF_EMAXDIMS
      INTEGER NF_ENAMEINUSE
      INTEGER NF_ENOTATT
      INTEGER NF_EMAXATTS
      INTEGER NF_EBADTYPE
      INTEGER NF_EBADDIM
      INTEGER NF_EUNLIMPOS
      INTEGER NF_EMAXVARS
      INTEGER NF_ENOTVAR
      INTEGER NF_EGLOBAL
      INTEGER NF_ENOTNC
      INTEGER NF_ESTS
      INTEGER NF_EMAXNAME
      INTEGER NF_EUNLIMIT
      INTEGER NF_ENORECVARS
      INTEGER NF_ECHAR
      INTEGER NF_EEDGE
      INTEGER NF_ESTRIDE
      INTEGER NF_EBADNAME
      INTEGER NF_ERANGE

      PARAMETER (NF_NOERR = 0)
      PARAMETER (NF_EBADID = -33)
      PARAMETER (NF_EEXIST = -35)
      PARAMETER (NF_EINVAL = -36)
      PARAMETER (NF_EPERM = -37)
      PARAMETER (NF_ENOTINDEFINE = -38)
      PARAMETER (NF_EINDEFINE = -39)
      PARAMETER (NF_EINVALCOORDS = -40)
      PARAMETER (NF_EMAXDIMS = -41)
      PARAMETER (NF_ENAMEINUSE = -42)
      PARAMETER (NF_ENOTATT = -43)
      PARAMETER (NF_EMAXATTS = -44)
      PARAMETER (NF_EBADTYPE = -45)
      PARAMETER (NF_EBADDIM = -46)
      PARAMETER (NF_EUNLIMPOS = -47)
      PARAMETER (NF_EMAXVARS = -48)
      PARAMETER (NF_ENOTVAR = -49)
      PARAMETER (NF_EGLOBAL = -50)
      PARAMETER (NF_ENOTNC = -51)
      PARAMETER (NF_ESTS = -52)
      PARAMETER (NF_EMAXNAME = -53)
      PARAMETER (NF_EUNLIMIT = -54)
      PARAMETER (NF_ENORECVARS = -55)
      PARAMETER (NF_ECHAR = -56)
      PARAMETER (NF_EEDGE = -57)
      PARAMETER (NF_ESTRIDE = -58)
      PARAMETER (NF_EBADNAME = -59)
      PARAMETER (NF_ERANGE = -60)

!
! Error handling modes:
!
      INTEGER NF_FATAL
      INTEGER NF_VERBOSE

      PARAMETER (NF_FATAL = 1)
      PARAMETER (NF_VERBOSE = 2)

!
! MISCELLANEOUS ROUTINES:
!
      CHARACTER*80   NF_INQ_LIBVERS

      CHARACTER*80   NF_STRERROR
!     		          (INTEGER	       NCERR)

      LOGICAL        NF_ISSYSERR
!			  (INTEGER             NCERR)

!
! CONTROL ROUTINES:
!

      INTEGER         NF_CREATE
!		          (CHARACTER*(*)       PATH,
!		           INTEGER             CMODE,
!		           INTEGER             NCID)

      INTEGER         NF_OPEN
!		          (CHARACTER*(*)       PATH,
!		           INTEGER             MODE,
!		           INTEGER             NCID)

      INTEGER         NF_SET_FILL
!		          (INTEGER             NCID,
!			   INTEGER             FILLMODE,
!			   INTEGER             OLD_MODE)

      INTEGER         NF_REDEF
!			  (INTEGER             NCID)

      INTEGER         NF_ENDDEF
!			  (INTEGER             NCID)

      INTEGER         NF_SYNC
!			  (INTEGER             NCID)

      INTEGER         NF_ABORT
!			  (INTEGER             NCID)

      INTEGER         NF_CLOSE
!			  (INTEGER             NCID)

      INTEGER         NF_DELETE
!			  (CHARACTER*(*)       NCID)

!
! GENERAL INQUIRY ROUTINES:
!

      INTEGER         NF_INQ
!			  (INTEGER             NCID,
!			   INTEGER             NDIMS,
!			   INTEGER             NVARS,
!			   INTEGER             NGATTS,
!			   INTEGER             UNLIMDIMID)

      INTEGER	      NF_INQ_NDIMS
!			  (INTEGER	       NCID,
!			   INTEGER	       NDIMS)

      INTEGER	      NF_INQ_NVARS
!			  (INTEGER             NCID,
!			   INTEGER             NVARS)

      INTEGER	      NF_INQ_NATTS
!			  (INTEGER             NCID,
!		           INTEGER             NGATTS)

      INTEGER         NF_INQ_UNLIMDIM
!			  (INTEGER             NCID,
!			   INTEGER             UNLIMDIMID)

!
! DIMENSION ROUTINES:
!

      INTEGER	      NF_DEF_DIM
!			  (INTEGER             NCID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             LEN,
!			   INTEGER             DIMID)

      INTEGER	      NF_INQ_DIMID
!			  (INTEGER             NCID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             DIMID)

      INTEGER         NF_INQ_DIM
!			  (INTEGER             NCID,
!			   INTEGER             DIMID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             LEN)

      INTEGER         NF_INQ_DIMNAME
!			  (INTEGER             NCID,
!			   INTEGER             DIMID,
!			   CHARACTER(*)        NAME)

      INTEGER	      NF_INQ_DIMLEN
!			  (INTEGER             NCID,
!			   INTEGER             DIMID,
!			   INTEGER             LEN)

      INTEGER         NF_RENAME_DIM
!			  (INTEGER             NCID,
!			   INTEGER             DIMID,
!			   CHARACTER(*)        NAME)

!
! GENERAL ATTRIBUTE ROUTINES:
!

      INTEGER         NF_INQ_ATT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             XTYPE,
!			   INTEGER             LEN)

      INTEGER	      NF_INQ_ATTID
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             ATTNUM)

      INTEGER         NF_INQ_ATTTYPE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             XTYPE)

      INTEGER	      NF_INQ_ATTLEN
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             LEN)

      INTEGER         NF_INQ_ATTNAME
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             ATTNUM,
!			   CHARACTER(*)        NAME)

      INTEGER         NF_COPY_ATT
!			  (INTEGER             NCID_IN,
!			   INTEGER             VARID_IN,
!			   CHARACTER(*)        NAME,
!			   INTEGER             NCID_OUT,
!			   INTEGER             VARID_OUT)

      INTEGER         NF_RENAME_ATT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        CURNAME,
!			   CHARACTER(*)        NEWNAME)

      INTEGER         NF_DEL_ATT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME)

!
! ATTRIBUTE PUT/GET ROUTINES:
!

      INTEGER         NF_PUT_ATT_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER	       LEN,
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_GET_ATT_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_PUT_ATT_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             XTYPE,
!			   INTEGER             LEN,
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_GET_ATT_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_PUT_ATT_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             XTYPE,
!			   INTEGER             LEN,
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_GET_ATT_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_PUT_ATT_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             XTYPE,
!			   INTEGER             LEN,
!			   INTEGER             IVALS(1))

      INTEGER         NF_GET_ATT_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             IVALS(1))

      INTEGER         NF_PUT_ATT_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             XTYPE,
!			   INTEGER             LEN,
!			   REAL                RVALS(1))

      INTEGER         NF_GET_ATT_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   REAL                RVALS(1))

      INTEGER         NF_PUT_ATT_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             XTYPE,
!			   INTEGER             LEN,
!			   DOUBLE              DVALS(1))

      INTEGER         NF_GET_ATT_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   DOUBLE              DVALS(1))

!
! GENERAL VARIABLE ROUTINES:
!

      INTEGER         NF_DEF_VAR
!			  (INTEGER             NCID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             DATATYPE,
!			   INTEGER             NDIMS,
!			   INTEGER             DIMIDS(1),
!			   INTEGER             VARID)

      INTEGER         NF_INQ_VAR
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             DATATYPE,
!			   INTEGER             NDIMS,
!			   INTEGER             DIMIDS(1),
!			   INTEGER             NATTS)

      INTEGER         NF_INQ_VARID
!			  (INTEGER             NCID,
!			   CHARACTER(*)        NAME,
!			   INTEGER             VARID)

      INTEGER         NF_INQ_VARNAME
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME)

      INTEGER         NF_INQ_VARTYPE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             XTYPE)

      INTEGER         NF_INQ_VARNDIMS
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             NDIMS)

      INTEGER         NF_INQ_VARDIMID
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             DIMIDS(1))

      INTEGER         NF_INQ_VARNATTS
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             NATTS)

      INTEGER         NF_RENAME_VAR
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        NAME)

      INTEGER         NF_COPY_VAR
!			  (INTEGER             NCID_IN,
!			   INTEGER             VARID,
!			   INTEGER             NCID_OUT)

!
! ENTIRE VARIABLE PUT/GET ROUTINES:
!

      INTEGER         NF_PUT_VAR_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_GET_VAR_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_PUT_VAR_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_GET_VAR_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_PUT_VAR_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_GET_VAR_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_PUT_VAR_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             IVALS(1))

      INTEGER         NF_GET_VAR_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             IVALS(1))

      INTEGER         NF_PUT_VAR_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   REAL                RVALS(1))

      INTEGER         NF_GET_VAR_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   REAL                RVALS(1))

      INTEGER         NF_PUT_VAR_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   DOUBLEPRECISION     DVALS(1))

      INTEGER         NF_GET_VAR_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   DOUBLEPRECISION     DVALS(1))

!
! SINGLE VARIABLE PUT/GET ROUTINES:
!

      INTEGER         NF_PUT_VAR1_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   CHARACTER*1         TEXT)

      INTEGER         NF_GET_VAR1_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   CHARACTER*1         TEXT)

      INTEGER         NF_PUT_VAR1_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   NF_INT1_T           I1VAL)

      INTEGER         NF_GET_VAR1_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   NF_INT1_T           I1VAL)

      INTEGER         NF_PUT_VAR1_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   NF_INT2_T           I2VAL)

      INTEGER         NF_GET_VAR1_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   NF_INT2_T           I2VAL)

      INTEGER         NF_PUT_VAR1_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   INTEGER             IVAL)

      INTEGER         NF_GET_VAR1_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   INTEGER             IVAL)

      INTEGER         NF_PUT_VAR1_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   REAL                RVAL)

      INTEGER         NF_GET_VAR1_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   REAL                RVAL)

      INTEGER         NF_PUT_VAR1_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   DOUBLEPRECISION     DVAL)

      INTEGER         NF_GET_VAR1_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             INDEX(1),
!			   DOUBLEPRECISION     DVAL)

!
! VARIABLE ARRAY PUT/GET ROUTINES:
!

      INTEGER         NF_PUT_VARA_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_GET_VARA_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_PUT_VARA_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_GET_VARA_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_PUT_VARA_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_GET_VARA_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_PUT_VARA_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER             IVALS(1))

      INTEGER         NF_GET_VARA_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER             IVALS(1))

      INTEGER         NF_PUT_VARA_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   REAL                RVALS(1))

      INTEGER         NF_GET_VARA_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   REAL                RVALS(1))

      INTEGER         NF_PUT_VARA_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   DOUBLEPRECISION     DVALS(1))

      INTEGER         NF_GET_VARA_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   DOUBLEPRECISION     DVALS(1))

!
! STRIDED VARIABLE PUT/GET ROUTINES:
!

      INTEGER         NF_PUT_VARS_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_GET_VARS_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_PUT_VARS_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_GET_VARS_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_PUT_VARS_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_GET_VARS_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_PUT_VARS_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER             IVALS(1))

      INTEGER         NF_GET_VARS_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER             IVALS(1))

      INTEGER         NF_PUT_VARS_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   REAL                RVALS(1))

      INTEGER         NF_GET_VARS_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   REAL                RVALS(1))

      INTEGER         NF_PUT_VARS_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   DOUBLEPRECISION     DVALS(1))

      INTEGER         NF_GET_VARS_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   DOUBLEPRECISION     DVALS(1))

!
! MAPPED VARIABLE PUT/GET ROUTINES:
!

      INTEGER         NF_PUT_VARM_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_GET_VARM_TEXT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   CHARACTER(*)        TEXT)

      INTEGER         NF_PUT_VARM_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_GET_VARM_INT1
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   NF_INT1_T           I1VALS(1))

      INTEGER         NF_PUT_VARM_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_GET_VARM_INT2
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   NF_INT2_T           I2VALS(1))

      INTEGER         NF_PUT_VARM_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   INTEGER             IVALS(1))

      INTEGER         NF_GET_VARM_INT
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   INTEGER             IVALS(1))

      INTEGER         NF_PUT_VARM_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   REAL                RVALS(1))

      INTEGER         NF_GET_VARM_REAL
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   REAL                RVALS(1))

      INTEGER         NF_PUT_VARM_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   DOUBLEPRECISION     DVALS(1))

      INTEGER         NF_GET_VARM_DOUBLE
!			  (INTEGER             NCID,
!			   INTEGER             VARID,
!			   INTEGER             START(1),
!			   INTEGER             COUNT(1),
!			   INTEGER	       STRIDE(1),
!			   INTEGER	       IMAP(1),
!			   DOUBLEPRECISION     DVALS(1))

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Begin netCDF 2.4 backward compatibility:
!

!      
! Functions in the FORTRAN interface
!
      INTEGER NCCRE
      INTEGER NCOPN
      INTEGER NCDDEF
      INTEGER NCDID
      INTEGER NCVDEF
      INTEGER NCVID
      INTEGER NCTLEN
      INTEGER NCSFIL


      INTEGER NCRDWR
      INTEGER NCCREAT
      INTEGER NCEXCL
      INTEGER NCINDEF
      INTEGER NCNSYNC
      INTEGER NCHSYNC
      INTEGER NCNDIRTY
      INTEGER NCHDIRTY
      INTEGER NCLINK
      INTEGER NCNOWRIT
      INTEGER NCWRITE
      INTEGER NCCLOB
      INTEGER NCNOCLOB
      INTEGER NCGLOBAL
      INTEGER NCFILL
      INTEGER NCNOFILL
      INTEGER MAXNCOP
      INTEGER MAXNCDIM
      INTEGER MAXNCATT
      INTEGER MAXNCVAR
      INTEGER MAXNCNAM
      INTEGER MAXVDIMS
      INTEGER NCNOERR
      INTEGER NCEBADID
      INTEGER NCENFILE
      INTEGER NCEEXIST
      INTEGER NCEINVAL
      INTEGER NCEPERM
      INTEGER NCENOTIN
      INTEGER NCEINDEF
      INTEGER NCECOORD
      INTEGER NCEMAXDS
      INTEGER NCENAME
      INTEGER NCENOATT
      INTEGER NCEMAXAT
      INTEGER NCEBADTY
      INTEGER NCEBADD
      INTEGER NCESTS
      INTEGER NCEUNLIM
      INTEGER NCEMAXVS
      INTEGER NCENOTVR
      INTEGER NCEGLOB
      INTEGER NCENOTNC
      INTEGER NCFOOBAR
      INTEGER NCSYSERR
      INTEGER NCFATAL
      INTEGER NCVERBOS
      INTEGER NCENTOOL

!
! netCDF data types:
!
      INTEGER NCBYTE
      INTEGER NCCHAR
      INTEGER NCSHORT
      INTEGER NCLONG
      INTEGER NCFLOAT
      INTEGER NCDOUBLE

      PARAMETER(NCBYTE = 1)
      PARAMETER(NCCHAR = 2)
      PARAMETER(NCSHORT = 3)
      PARAMETER(NCLONG = 4)
      PARAMETER(NCFLOAT = 5)
      PARAMETER(NCDOUBLE = 6)

!     
!     masks for the struct NC flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      PARAMETER(NCRDWR = 1)
!     in create phase, cleared by ncendef 
      PARAMETER(NCCREAT = 2)
!     on create destroy existing file 
      PARAMETER(NCEXCL = 4)
!     in define mode, cleared by ncendef 
      PARAMETER(NCINDEF = 8)
!     synchronise numrecs on change (X'10')
      PARAMETER(NCNSYNC = 16)
!     synchronise whole header on change (X'20')
      PARAMETER(NCHSYNC = 32)
!     numrecs has changed (X'40')
      PARAMETER(NCNDIRTY = 64)	
!     header info has changed (X'80')
      PARAMETER(NCHDIRTY = 128)
!     prefill vars on endef and increase of record, the default behavior
      PARAMETER(NCFILL = 0)
!     do not fill vars on endef and increase of record (X'100')
      PARAMETER(NCNOFILL = 256)
!     isa link (X'8000')
      PARAMETER(NCLINK = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      PARAMETER(NCNOWRIT = 0)
      PARAMETER(NCWRITE = NCRDWR)
      PARAMETER(NCCLOB = NF_CLOBBER)
      PARAMETER(NCNOCLOB = NF_NOCLOBBER)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      INTEGER NCUNLIM
      PARAMETER(NCUNLIM = 0)

!     
!     attribute id to put/get a global attribute
!     
      PARAMETER(NCGLOBAL  = 0)

!     
!     Advisory Maximums:
!     
      PARAMETER(MAXNCOP = 32)
      PARAMETER(MAXNCDIM = 100)
      PARAMETER(MAXNCATT = 2000)
      PARAMETER(MAXNCVAR = 2000)
!     Not enforced 
      PARAMETER(MAXNCNAM = 128)
      PARAMETER(MAXVDIMS = MAXNCDIM)

!     
!     Global netcdf error status variable
!     Initialized in error.c
!     

!     No Error 
      PARAMETER(NCNOERR = 0)
!     Not a netcdf id 
      PARAMETER(NCEBADID = 1)
!     Too many netcdfs open 
      PARAMETER(NCENFILE = 2)	
!     netcdf file exists && NCNOCLOB
      PARAMETER(NCEEXIST = 3)
!     Invalid Argument 
      PARAMETER(NCEINVAL = 4)
!     Write to read only 
      PARAMETER(NCEPERM = 5)
!     Operation not allowed in data mode 
      PARAMETER(NCENOTIN = 6)	
!     Operation not allowed in define mode 
      PARAMETER(NCEINDEF = 7)	
!     Coordinates out of Domain 
      PARAMETER(NCECOORD = 8)
!     MAXNCDIMS exceeded 
      PARAMETER(NCEMAXDS = 9)
!     String match to name in use 
      PARAMETER(NCENAME = 10)	
!     Attribute not found 
      PARAMETER(NCENOATT = 11)
!     MAXNCATTRS exceeded 
      PARAMETER(NCEMAXAT = 12)
!     Not a netcdf data type 
      PARAMETER(NCEBADTY = 13)
!     Invalid dimension id 
      PARAMETER(NCEBADD = 14)	
!     NCUNLIMITED in the wrong index 
      PARAMETER(NCEUNLIM = 15)
!     MAXNCVARS exceeded 
      PARAMETER(NCEMAXVS = 16)
!     Variable not found 
      PARAMETER(NCENOTVR = 17)
!     Action prohibited on NCGLOBAL varid 
      PARAMETER(NCEGLOB = 18)
!     Not a netcdf file 
      PARAMETER(NCENOTNC = 19)
      PARAMETER(NCESTS = 20)
      PARAMETER (NCENTOOL = 21)	
      PARAMETER(NCFOOBAR = 32)
      PARAMETER(NCSYSERR = -1)

!     
!     Global options variable. Used to determine behavior of error handler.
!     Initialized in lerror.c
!     
      PARAMETER(NCFATAL = 1)
      PARAMETER(NCVERBOS = 2)

!
!     Default fill values.  These must be the same as in the C interface.
!
      INTEGER FILBYTE
      INTEGER FILCHAR
      INTEGER FILSHORT
      INTEGER FILLONG
      REAL FILFLOAT
      DOUBLEPRECISION FILDOUB

      PARAMETER (FILBYTE = -127)
      PARAMETER (FILCHAR = 0)
      PARAMETER (FILSHORT = -32767)
      PARAMETER (FILLONG = -2147483647)
      PARAMETER (FILFLOAT = 9.9692099683868690e+36)
      PARAMETER (FILDOUB = 9.9692099683868690e+36)
