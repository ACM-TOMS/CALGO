// Time-stamp:  <2015-12-18 14:57:58 m>
#include "csample.h"
void SAMPLE_(struct SAMPLE_TY_ *s)
{
  int idat[2];
  CTYPE_ rdat[1] = {(fp)3.1415926535897932384626433832795028841971693993};
  
  switch (s->what) {
  case setup_sample: // The setup call
    strcpy(s->se.ename, SAMPLE_NAME_); 
    return;
  case partial_message: // Partial print of error1
    idat[0]=1; idat[1]=s->what;
    MESSY_(&s->se, "$E34Start of error message from $Psample what=$I.$N"
           "I want some $R!$C", m_idat,2,idat, m_rdat,1,rdat,m_ptext,qplet_, 0);
    return;
  case finish_message: // Finish this error message
    idat[0] = s->what;
    MESSY_(&s->se, "Finish with what=$I.", m_idat,1,idat, 0);
    return;
  default:
    idat[0] = 99; idat[1]= s->what;
    MESSY_(&s->se, "$E56Fatal error called " SAMPLE_NAME_ " with what = $I.",
           m_idat,2,idat, 0);
    return;
  }
}
