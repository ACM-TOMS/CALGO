package com.googlecode.compressingcircularbuffers;
/**
 * A simple version of a compressing buffer. Not suitable for actual use
 * (lacks many required methods). Used only for comparative performance
 * testing with CompressingCircularBuffer. 
 * 
 */
class SimpleCompressingBuffer {
    private double[] buffer = null; 
    private CompressionMethod compressionMethod = null; 
    private int iElement = -1;  
    private int nElements = 0;  
    private int compressionMultiplier = 1; 
    private int minCompressionRatio = 1;
    private int maxCompressionRatio = Integer.MAX_VALUE;
    
    public SimpleCompressingBuffer(int nElements, CompressionMethod compressionMethod, 
                                     int minCompressionRatio, int maxCompressionRatio) {
      if (nElements < 1 || nElements % 2 == 0)
        throw new IllegalArgumentException("nElements=" + nElements +
                                           "; nElements must be odd and at least 1.");
      else if (minCompressionRatio < 1)
        throw new IllegalArgumentException("minCompressionRatio (" +
                                           minCompressionRatio + ") must be at least 1 ");
      else if (maxCompressionRatio < minCompressionRatio)
        throw new IllegalArgumentException("maxCompressionRatio (" +
                                           maxCompressionRatio + ") must be >= " +
                                          " minCompressionRatio (" + minCompressionRatio +")");
      this.compressionMethod = compressionMethod;
      this.nElements = nElements;
      this.minCompressionRatio = minCompressionRatio;
      this.maxCompressionRatio = maxCompressionRatio;
      buffer = new double[nElements];
   }
    
   public SimpleCompressingBuffer(int nElements, CompressionMethod compressionMethod,
                                     int minCompressionRatio) {
      this(nElements, compressionMethod, minCompressionRatio, Integer.MAX_VALUE);
   }
   public SimpleCompressingBuffer(int nElements, CompressionMethod compressionMethod) {
     this(nElements, compressionMethod, 1, Integer.MAX_VALUE);
   }
   public SimpleCompressingBuffer(int nElements) {
     this(nElements, new AveragingCompressionMethod(),  1, Integer.MAX_VALUE);
   }
   public SimpleCompressingBuffer() {
     this(1, new AveragingCompressionMethod(),  1, Integer.MAX_VALUE);
   }

   public void update(double sample) {
      compressionMethod.update(sample);
      if (compressionMethod.getCount() == compressionMultiplier*minCompressionRatio) {
        // next compressed value is ready to be stored
        double newData = compressionMethod.getCompressedValue();
        iElement++;
        if (iElement >= nElements &&
            compressionMultiplier * minCompressionRatio > maxCompressionRatio/2) {
          // max compression ratio exceeded; treat like a circular buffer
          buffer[iElement % nElements] = newData;
        }
        else if (iElement == nElements) { // buffer is full at this compression ratio, so...
          compressionMultiplier *= 2;     // double ratio to double the capacity
          for (int i = 1; i < nElements-1; i+=2) 
            buffer[i/2] = compressionMethod.mergeValues(buffer[i-1], buffer[i]);
          // since capacity is doubled, we are now on the middle element: 
          iElement = (nElements-1)/2;
          double older = buffer[nElements-1];
          buffer[iElement] = compressionMethod.mergeValues(older, newData);;
        }
        else {
          buffer[iElement] = newData;
        }
        compressionMethod.reset(); 
      }
   }

   /**
    * Returns the simple circular buffer to it's initial state,
    * discarding all previously seen data.
    * 
    */
   public void reset() {
     iElement = -1;  
     compressionMultiplier = 1;
     compressionMethod.reset();
   }
   
}
