
#ifndef __sdpa_block_h__
#define __sdpa_block_h__

#include "sdpa_include.h"

namespace sdpa {

class BlockStruct
{
public:
  enum BlockType {btSDP,btSOCP,btLP};
  int  nBlock;
  int* blockStruct;
  int* blockNumber;
  BlockType* blockType;
  int  SDP_nBlock;
  int* SDP_blockStruct;
  int  SOCP_nBlock;
  int* SOCP_blockStruct;
  int  LP_nBlock;

  BlockStruct();
  ~BlockStruct();
  void initialize(int nBlock);
  void terminate();
  void makeInternalStructure();
  void display(FILE* fpOut = stdout);

};

} //end of namespace 'sdpa'

#endif // __sdpa_block_h__
