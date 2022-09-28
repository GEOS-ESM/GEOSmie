      INTEGER*4, PARAMETER:: NANG_MAX=1000, NG_MAX=100000
C**********************************************************************
      REAL*8, PARAMETER  :: PI=DACOS(-1d0)
      REAL*8, PARAMETER  :: D2R=PI/180D0
      REAL*8, PARAMETER  :: R2D=180D0/PI
C**********************************************************************
      INTEGER, PARAMETER :: MAXABS = 1, MAXPERCENT = 2, MSRE = 3
      INTEGER, PARAMETER :: ERRTYP =  MAXABS ! MAXABS, MAXPERCENT or MSRE
      REAL*8, PARAMETER  :: ang_min = 0.0D0*D2R
      REAL*8, PARAMETER  :: ang_max = 180.0D0*D2R
      REAL*8, PARAMETER  :: DESIRED_ERR = 0.02D0
      INTEGER, PARAMETER :: DELTA_NG = 10
      INTEGER, PARAMETER :: USE_ALT_ANG = 1
      INTEGER, PARAMETER :: NSPHER = 129
      INTEGER, PARAMETER :: MIN_NSPHER = 129
      

