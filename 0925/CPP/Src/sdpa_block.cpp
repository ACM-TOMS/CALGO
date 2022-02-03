
#include "sdpa_block.h"

namespace sdpa {
  
BlockStruct::BlockStruct()
{
  nBlock      = 0;
  blockStruct = NULL;
  blockNumber = NULL;
  blockType   = NULL;
  SDP_nBlock       = 0;
  SDP_blockStruct  = NULL;
  SOCP_nBlock      = 0;
  SOCP_blockStruct = NULL;
  LP_nBlock        = 0;
}

BlockStruct::~BlockStruct()
{
  terminate();
}

void BlockStruct::initialize(int nBlock)
{
  this->nBlock = nBlock;
  NewArray(blockStruct,int,nBlock);
  NewArray(blockType, BlockType, nBlock);
  NewArray(blockNumber,int,nBlock);
  SDP_nBlock       = 0;
  SDP_blockStruct  = NULL;
  SOCP_nBlock      = 0;
  SOCP_blockStruct = NULL;
  LP_nBlock        = 0;
}

void BlockStruct::terminate()
{
  DeleteArray(blockStruct);
  DeleteArray(blockNumber);
  DeleteArray(blockType);
  DeleteArray(SDP_blockStruct);
  DeleteArray(SOCP_blockStruct);
}

void BlockStruct::makeInternalStructure()
{
  SDP_nBlock  = 0;
  SOCP_nBlock = 0;
  LP_nBlock   = 0;
  for (int l=0; l<nBlock; l++){
    #if 0
    rMessage("blockStruct[" << l << "] = "<<  blockStruct[l]
	     << ": blockType[" << l << "] = " << blockType[l]);
    #endif
    if (blockStruct[l] >= 2 && blockType[l] == btSDP) {
      blockType[l]   = btSDP;
      blockNumber[l] = SDP_nBlock;
      SDP_nBlock++;
    } else if (blockStruct[l] < 0 || blockType[l] == btLP) {
      blockType[l]   = btLP;
      if (blockStruct[l] < 0) {
	blockStruct[l] = - blockStruct[l];
      }
      blockNumber[l] = LP_nBlock;
      LP_nBlock += blockStruct[l];
    } else if (blockStruct[l] == 1) {
      blockType[l]   = btLP;
      blockStruct[l] = 1;
      blockNumber[l] = LP_nBlock;
      LP_nBlock += blockStruct[l];
    } else {
      rError("block struct");
    }
  }
  

  if (SDP_nBlock > 0) {
    NewArray(SDP_blockStruct, int,SDP_nBlock);
  }
  SDP_nBlock = 0;
  for (int l=0; l<nBlock; l++){
    if (blockType[l] == btSDP) {
      SDP_blockStruct[SDP_nBlock] = blockStruct[l];
      SDP_nBlock++;
    } 
  }

  #if 0
  if (SOCP_nBlock > 0) {
    NewArray(SOCP_blockStruct,int,SOCP_nBlock);
  }
  SOCP_nBlock = 0;
  for (int l=0; l<nBlock; l++){
    if (blockType[l] == btSOCP) {
      SOC_blockStruct[SOCP_nBlock] = blockStruct[l];
      SOCP_nBlock++;
    } 
  }
  #endif
  
}

void BlockStruct::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"--- BlockStruct ---\n");
  fprintf(fpout,"nBlock = %d\n",nBlock);
  fprintf(fpout,"blockStruct = \n");
  for (int l=0; l<nBlock; ++l) {
    fprintf(fpout,"%5d,",blockStruct[l]);
  }
  fprintf(fpout,"\n");
  fprintf(fpout,"blockNumber = \n");
  for (int l=0; l<nBlock; ++l) {
    fprintf(fpout,"%5d,",blockNumber[l]);
  }
  fprintf(fpout,"\n");
  fprintf(fpout,"blockType = \n");
  for (int l=0; l<nBlock; ++l) {
    char displaychar = '-';
    if (blockType[l] == btSDP) {
      displaychar= 'S';
    } else if (blockType[l] == btSOCP) {
      displaychar= 'Q';
    } else if (blockType[l] == btLP) {
      displaychar= 'L';
    }
    fprintf(fpout,"    %c,",displaychar);
  }
  fprintf(fpout,"\n");
  fprintf(fpout,"SDP_nBlock = %d\n",SDP_nBlock);
  fprintf(fpout,"SDP_blockStruct = \n");
  for (int l=0; l<SDP_nBlock; ++l) {
    fprintf(fpout,"%5d,",SDP_blockStruct[l]);
  }
  fprintf(fpout,"\n");
  fprintf(fpout,"SOCP_nBlock = %d\n",SOCP_nBlock);
  fprintf(fpout,"SOCP_blockStruct = \n");
  for (int l=0; l<SOCP_nBlock; ++l) {
    fprintf(fpout,"%5d,",SOCP_blockStruct[l]);
  }
  fprintf(fpout,"\n");
  fprintf(fpout,"LP_nBlock = %d\n",LP_nBlock);
  fprintf(fpout,"--- BlockStruct ---\n");
}


} // end of namespace 'sdpa'
