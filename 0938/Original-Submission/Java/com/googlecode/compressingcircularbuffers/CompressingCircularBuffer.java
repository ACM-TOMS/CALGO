package com.googlecode.compressingcircularbuffers;
/**
 * Maintains a compressed, fixed-memory-footprint, "sketch" (summary) of
 * a sampled data sequence. The compressed sequence is formed by
 * applying a specified compression method to a series of approximately
 * equal sized contiguous sequential chunks of the original sequence.
 * <p>
 *
 * For example, with the default, <tt>AveragingCompressionMethod</tt>,
 * sequential raw data samples (fed in via the <tt>update</tt> method)
 * are averaged together in approximately equal sized contiguous chunks
 * to form the compressed sequence. As new data comes in, those chunks
 * expand in size (with associated averages adjusted accordingly) so
 * that a buffer of the specified size can still contain a compressed
 * form of the growing sequence.  <p>
 *
 * In general, these chunk sizes will either <tt>2<sup>k</sup></tt> or
 * <tt>2<sup>k+1</sup></tt>, where <tt>k</tt> is chosen so that the aggregated
 * samples fill up the available buffer space (<tt>k</tt> starts at zero
 * and increases as more samples come in and the buffer becomes full at
 * the current compression level). Chunks at the beginning or at the end of
 * the sequence use the larger aggregation interval, with the smaller
 * interval used in the middle. You can use the methods
 * <tt>getFirstSampleIndex</tt> and <tt>getLastSampleIndex</tt> to
 * determine the number of samples aggregated into any specific
 * compressed value.  <p>
 * 
 * Performance is very similar to that of an ordinary circular buffer:
 * The algorithm updates the compressed sequence in <i>at most</i> one
 * read from, and two writes to, its internal buffer per update.  Thus
 * it is appropriate for use when update times must be very predictable,
 * since there is never any 'stalling' as large chunks of memory are
 * reorganized to reflect a change in compression level.  For example,
 * if 'memory' were a large file on a hard disk, the algorithm could
 * still have acceptable worst case update times (a couple of seeks/writes
 * per update). As with ordinary circular buffers, it is also close to
 * optimal in memory use: all buffer elements are filled with useful
 * data.  <p>
 * 
 * This class uses an indexing scheme that exploits basic properties of
 * modular arithmetic to quickly locate exactly those positions
 * containing the 'compressed away' data that are available to store the
 * new data. That's why it's more efficient than more straightforward
 * approaches. In the associated paper "Compressing Circular Buffers" I
 * explain this indexing scheme more completely.  <p>
 * 
 * Copyright, 2010, John C. Gunther. All rights reserved.  Released under the
 * terms of the Apache 2.0 licence. 
 * 
 */
public class CompressingCircularBuffer {

    private double[] buffer = null; // stores compressed version of sampled sequence
    private CompressionMethod compressionMethod = null; // thinning, averaging, etc.
    private int iElement = -1;  // index used to access last compressed value via:
                                // buffer[(iElement*compressionRatio) % nElements]
    private int nElements = 0;  // number of elements in the compressed sequence                     
// compressionMultiplier*minCompressionRatio = # points contained in new compressed values
    private int compressionMultiplier = 1; 
    private int minCompressionRatio = 1;
    /**
     * Upper bound on maximum number of updates that can be performed
     * without a reset. This bound arises because integer related indexing
     * calculations must remain in the valid range for integers.
     *
     */
    public final static int MAX_ALLOWED_UPDATES = Integer.MAX_VALUE/2;

    /**
     * Compressed element index returned when no compressed element
     * matching given criteria is available.
     *
     */
    public final static int NOSUCHINDEX = -Integer.MAX_VALUE;
    /**
     * Default upper bound on compression ratio if none is specified in
     * constructor.
     */ 
    public final static int DEFAULT_MAX_COMPRESSION_RATIO = Integer.MAX_VALUE/4;
    private int maxCompressionRatio = DEFAULT_MAX_COMPRESSION_RATIO;
    /**
     * Creates a new compressing circular buffer.
     * <p>
     * 
     * @param nElements number of compressed elements used to represent
     * sequence. Must be an odd integer and at least 1. Default is 1.
     *
     * @param compressionMethod how compression will be performed. Use
     * an instance of one of the predefined compression methods (<tt>ThinningCompressionMethod,
     * AveragingCompressionMethod, SummingCompressionMethod, MinimizingCompresionMethod,
     * MaximizingCompressionMethod, or RandomSamplingCompressionMethod</tt>) or else
     * create your own by implementing the <tt>CompressionMethod</tt> interface.
     * Default is <tt>new AveragingCompressionMethod()</tt>.
     *
     * @param minCompressionRatio minimum number of data points
     * contained in any compressed value. Default is 1.
     * 
     * @param maxCompressionRatio maximum number of data points in any
     * compressed value. Must be <tt>&ge; minCompressionRatio</tt>.  Note
     * that compressed values always contain a number of points equal to
     * <tt>minCompressionRatio*2<sup>k</sup></tt> with <tt>k</tt> a non-negative integer.
     * Thus actual maximum is the largest integer of this form that does
     * not exceed the <tt>maxCompressionRatio</tt> that you specify. Once this
     * limit is reached, the sequence will behave like an ordinary
     * circular buffer (wrapping around and overwriting the oldest
     * compressed value). Default is <tt>DEFAULT_MAX_COMPRESSION_RATIO</tt>.
     *
     * @see #DEFAULT_MAX_COMPRESSION_RATIO DEFAULT_MAX_COMPRESSION_RATIO
     * @see ThinningCompressionMethod ThinningCompressionMethod
     * @see AveragingCompressionMethod AveragingCompressionMethod
     * @see SummingCompressionMethod SummingCompressionMethod
     * @see MinimizingCompressionMethod MinimizingCompressionMethod
     * @see MaximizingCompressionMethod MaximizingCompressionMethod
     * @see RandomSamplingCompressionMethod RandomSamplingCompressionMethod
     * 
     */
    public CompressingCircularBuffer(int nElements, CompressionMethod compressionMethod, 
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
      this.nElements = nElements;
      this.compressionMethod = compressionMethod;
      this.minCompressionRatio = minCompressionRatio;
      this.maxCompressionRatio = maxCompressionRatio;
      buffer = new double[nElements];
   }
    
    /**
     * Creates a new compressing circular buffer.
     * <p>
     * 
     * Convenience constructor. Same as:
     *
     * <pre> 
     * CompressingCircularBuffer(nElements, compressionMethod,
     *                           minCompressionRatio, DEFAULT_MAX_COMPRESSION_RATIO)
     * </pre>
     *
     * @see #CompressingCircularBuffer(int, CompressionMethod, int, int) 
     *   CompressingCircularBuffer(int, CompressionMethod, int, int)
     *   
     */
   public CompressingCircularBuffer(int nElements, CompressionMethod compressionMethod,
                                     int minCompressionRatio) {
      this(nElements, compressionMethod, minCompressionRatio, DEFAULT_MAX_COMPRESSION_RATIO);
   }
    /**
     * Creates a new compressing circular buffer.
     * <p>
     * 
     * Convenience constructor. Same as:
     *
     * <pre> 
     * CompressingCircularBuffer(nElements, compressionMethod,
     *                           1, DEFAULT_MAX_COMPRESSION_RATIO)
     * </pre>
     * 
     * @see #CompressingCircularBuffer(int, CompressionMethod, int, int) 
     *   CompressingCircularBuffer(int, CompressionMethod, int, int)
     *
     */
   public CompressingCircularBuffer(int nElements, CompressionMethod compressionMethod) {
     this(nElements, compressionMethod, 1, DEFAULT_MAX_COMPRESSION_RATIO);
   }
    /**
     * Creates a new compressing circular buffer.
     * <p>
     * 
     * Convenience constructor. Same as:
     *
     * <pre> 
     * CompressingCircularBuffer(nElements, new AveragingCompressionMethod(),
     *                           1, DEFAULT_MAX_COMPRESSION_RATIO)
     * </pre>
     * 
     * @see #CompressingCircularBuffer(int, CompressionMethod, int, int) 
     *   CompressingCircularBuffer(int, CompressionMethod, int, int)
     *
     */
   public CompressingCircularBuffer(int nElements) {
     this(nElements, new AveragingCompressionMethod(),  1, DEFAULT_MAX_COMPRESSION_RATIO);
   }
    /**
     * Creates a new compressing circular buffer.
     * <p>
     * 
     * Convenience constructor. Same as:
     *
     * <pre> 
     * CompressingCircularBuffer(1, new AveragingCompressionMethod(),
     *                           1, DEFAULT_MAX_COMPRESSION_RATIO)
     * </pre>
     * 
     * @see #CompressingCircularBuffer(int, CompressionMethod, int, int) 
     *   CompressingCircularBuffer(int, CompressionMethod, int, int)
     *
     */
   public CompressingCircularBuffer() {
     this(1, new AveragingCompressionMethod(),  1, DEFAULT_MAX_COMPRESSION_RATIO);
   }
/**
 * Updates the compressed representation of the sequence to include a newly
 * collected data value (sequence element). 
 * <p>
 * 
 * For example, if you were sampling an on-line sensor once per minute,
 * you might call update every minute, feeding in each newly sampled
 * value. In that case, the compressing circular buffer would maintain
 * an up-to-date compressed "sketch" of the entire sequence; that sketch
 * would fit into the number of buffer elements that you specified in
 * the compressing circular buffer's constructor.  <p>
 * 
 * @param sample next value of a sampled sequence to be incorporated into the compressed representation.
 * 
 */
   public void update(double sample) {
      compressionMethod.update(sample);
      if (compressionMethod.getCount() == compressionMultiplier*minCompressionRatio) {
      // current compressed value has enough data in it to be stored
        double newData = compressionMethod.getCompressedValue();
        iElement++;
        if (iElement * compressionMultiplier * minCompressionRatio > MAX_ALLOWED_UPDATES)
          throw new IllegalStateException(
"Number of updates has exceeded MAX_ALLOWED_UPDATES (" + MAX_ALLOWED_UPDATES + "). " +
"This maximum is enforced to assure that Integer.MAX_VALUE is never be exceeded; such an " + 
"integer overflow would make CompressingCircularBuffer's index-related calculations invalid.");                                         
        else if (iElement == nElements &&
            compressionMultiplier * minCompressionRatio  <= maxCompressionRatio/2) {
          // buffer is full, and we are not at max compression ratio
          compressionMultiplier *= 2;     // double ratio to double the capacity
          // since capacity is doubled, we are now on the middle element
          // (it was the final element at the previous compression level)
          iElement = (nElements-1)/2;
          // the first 'half' of the middle element is already in the buffer;
          // the second 'half' is the newly compressed data value
          int index = (iElement*compressionMultiplier) % nElements;
          double older = buffer[index];
          double mergedValue = compressionMethod.mergeValues(older, newData);
          buffer[index] = mergedValue;
        }
        else if (iElement >= nElements) { // max compression exceeded, just wrap around
          int index = iElement*compressionMultiplier;
          index %= nElements;
          buffer[index] = newData;
        }
        else if (1 == compressionMultiplier) { // startup is special (no merging/compressing)
          buffer[iElement] = newData;
        }
        else { // normal case: merge in to-be-overwritten data; overwrite with new data
          int newerIndex = iElement * compressionMultiplier;
          int olderIndex = newerIndex - compressionMultiplier/2;
          newerIndex %= nElements;
          olderIndex %= nElements;
          double older = buffer[olderIndex];
          double newer = buffer[newerIndex];
          double merged = compressionMethod.mergeValues(older, newer);    
          buffer[olderIndex] = merged;
          buffer[newerIndex] = newData;
        }
        compressionMethod.reset(); 
      }
   }

   /**
    * Retrieves the index of the compressed value whose aggregation
    * interval includes the given raw data sample. 
    *
    * @param iRawDataSample raw data sample index
    * @return index of compressed value whose aggregation interval
    * includes the index of the given data sample, or
    * <tt>NOSUCHINDEX</tt> if no compressed value includes the given sample.
    * 
    *
    */
   public int getCompressedIndex(int iRawDataSample) {
     int result = NOSUCHINDEX;
     int iSample = iRawDataSample/minCompressionRatio;
     int iChunk = iSample/compressionMultiplier;

     if (iRawDataSample < 0)
       throw new IllegalArgumentException("iRawDataSample ("+iRawDataSample+") is < 0");
     else if (iRawDataSample >= getNCompressedSamples())
       return result;
     else if (iElement >= nElements) { 
// exceeded max compression level; acts much like an ordinary circular buffer
       if (iElement - iChunk < nElements) 
         result = nElements-1 - (iElement-iChunk); 
     }
     else if (1 == compressionMultiplier) { // startup, buffer may not be full
       if (iSample <= iElement) {
         result = iSample;
       }
     }
     else { // normal case.
       int midElement = (nElements - 1)/2;
       // Number of big 'hops' (raw data index increments between
       // successive compressed values) at the beginning of the
       // compressed data sequence. In general, hops at the beginning
       // and the end of the compressed sequence will be 'big hops' (that
       // is, hops spaced at the current compression ratio) whereas
       // compressed values in the middle are spaced at the previous
       // compression ratio (which is half the current ratio).
       int leadingBigHops = iElement - midElement;
       if (iChunk < leadingBigHops) // in full compressed, initial range
         result = iChunk;
       else if (iChunk < midElement) // in half-compressed, middle range
         // the twice as big hops at beginning get 'counted twice'
         result = iSample/(compressionMultiplier/2) - leadingBigHops;
       else { // in full compressed, final range
         int trailingBigHops = iChunk - midElement;
         // the twice as big hops at both ends get 'counted twice'
         result = 2*iChunk - leadingBigHops - trailingBigHops;
       }
     }
     return result;
   }
   /*
    * Retrieves the index of the first uncompressed sequence value included in the
    * compressed value with the specified index, expressed in units of
    * <tt>minCompressionRatio</tt>.
    * <p>
    *
    * Returns <tt>NOSUCHINDEX</tt> if there is currently no compressed value at the given
    * index. Does not check that iCompressed is in the required range;
    * that's done by the public wrapper methods.
    *
    *
    */
   private int getFirstSampleIndexUnchecked(int iCompressed) {
     int result = NOSUCHINDEX;
     
     if (iElement >= nElements) { 
// exceeded max compression level; acts much like an ordinary circular buffer
       result = iElement - ((nElements-1) - iCompressed);       
       result *= compressionMultiplier;
     }
     else if (1 == compressionMultiplier) { // startup, buffer may not be full
       if (iCompressed <= iElement) {
         result = iCompressed;
       }
     }
     else { // normal case.
       int midElement = (nElements - 1)/2;
       // Number of big 'hops' (raw data index increments between
       // successive compressed values) at the beginning of the
       // compressed data sequence. In general, hops at the beginning
       // and the end of the compressed sequence will be 'big hops' (that
       // is, hops spaced at the current compression ratio) whereas
       // compressed values in the middle are spaced at the previous
       // compression ratio (which is half the current ratio).
       int leadingBigHops = iElement - midElement;
       if (iCompressed <= leadingBigHops) // in full compressed, initial range
         result = iCompressed * compressionMultiplier;
       else if (iCompressed < nElements-leadingBigHops) // in half-compressed, middle range
         // the twice as big hops at beginning get 'counted twice'
         result = (iCompressed+leadingBigHops)*compressionMultiplier/2;
       else { // in full compressed, final range
         int trailingBigHops = iCompressed - midElement;
         // the twice as big hops at both ends get 'counted twice'
         result = (iCompressed + leadingBigHops + trailingBigHops)*compressionMultiplier/2;
       }
     }
     return result;
   }
   
   /**
    * Retrieves the index of the first uncompressed sequence value included in the
    * specified compressed value; uncompressed sequence indexes start at 0.
    * <p>
    *
    * Returns <tt>NOSUCHINDEX</tt> if there is currently no compressed value at the given
    * index.
    * <p>
    * 
    * This method can be used when you need to determine the range of underlying
    * uncompressed sequence values that a compressed value incorporates.
    * <p>
    * 
    * @param iCompressed the compressed sequence element whose sample index is to be retrieved
    * @return index of the first uncompressed sequence value incorporated into the specified compressed value.
    * 
    * @see #getLastSampleIndex(int) getLastSampleIndex
    * @see #NOSUCHINDEX NOSUCHINDEX
    */
   public int getFirstSampleIndex(int iCompressed) {
     if (iCompressed < 0 || iCompressed >= nElements)
       throw new IllegalArgumentException("iCompressed=" + iCompressed +
          "; iCompressed must be >= 0 and < nElements (which is " + nElements + ")");
     int result = getFirstSampleIndexUnchecked(iCompressed);
     if (result != NOSUCHINDEX)
       result *= minCompressionRatio;
     return result;
   }
   /**
    * Retrieves the index of the last uncompressed sequence value incorporated into the
    * specified compressed value.
    * <p>
    *
    * Returns <tt>NOSUCHINDEX</tt> if there is currently no compressed value at the given
    * index.
    *
    * @param iCompressed the compressed sequence element whose sample index is to be retrieved
    * @return index of the last uncompressed sequence value incorporated into the
    * specified compressed value.
    * 
    * @see #getFirstSampleIndex(int) getFirstSampleIndex
    * @see #NOSUCHINDEX NOSUCHINDEX
    * 
    */
   public int getLastSampleIndex(int iCompressed) {
     if (iCompressed < 0 || iCompressed >= nElements)
       throw new IllegalArgumentException("iCompressed=" + iCompressed +
                                          "; iCompressed must be >= 0 and < nElements (which is " + nElements + ")");
     int result;
     if (1 == compressionMultiplier)
       result = getFirstSampleIndex(iCompressed);
     else if (iCompressed == nElements-1)
       result = getFirstSampleIndexUnchecked(iCompressed) + compressionMultiplier - 1;
     else 
       result = getFirstSampleIndexUnchecked(iCompressed+1)-1;
     
     if (result != NOSUCHINDEX)
       result *= minCompressionRatio;
     
     return result;
   }

   /**
    *
    * Retrieves the number of samples fed to the circular buffer that
    * have been incorporated into the compressed sequence. The most
    * recently collected samples that have not yet been incorporated
    * into the compressed sequence are not included in the count.
    * 
    * @return number of samples fed into the compressing circular buffer since
    *   it was created (or since the last <tt>reset</tt>), and that have
    *   been incorporated into the compressed sequence.
    */
   public int getNCompressedSamples() {
     int result;
     if (compressionMultiplier == 1)
       result = minCompressionRatio * (iElement + 1);
     else
       result = minCompressionRatio *
                (getFirstSampleIndexUnchecked(nElements-1) + compressionMultiplier) ;
     return result;
   }

   /**
    *
    * Retrieves the number of samples fed to the circular buffer so far.
    * In other words, the number of times that update was called since
    * instantiation or reset.
    * 
    * @return number of samples fed into the compressing circular buffer since
    *   it was created (or since the last <tt>reset</tt>).
    */
   public int getNSamples() {
     int result = getNCompressedSamples();
     result += compressionMethod.getCount();
     return result;
   }
   /**
    *
    * Retrieves the number of compressed values within the circular
    * buffer. 
    * <p>
    * 
    * <i>Note:</i> Except during startup, this will always equal the number
    * of elements in the buffer specified in the constructor.
    * 
    * @return number of compressed elements in the buffer
    */
   public int getNCompressed() {
     int result;
     if (compressionMultiplier == 1)
       result = iElement + 1;
     else
       result = nElements;

     return result;
   }

   /**
    * Retrieves the number of points that will be included in the next
    * new compressed value that will be formed from samples fed in
    * via the <tt>update</tt> method.
    * <p>
    *
    * Initially, this starts at <tt>minCompressionRatio</tt>, until the
    * compressed values buffer fills up, in which case it doubles, until
    * the buffer is filled again at the higher compression ratio, then
    * it doubles again, and so forth.
    *
    * @return number of samples (updates) included in each newly added compressed 
    *   sequence value.
    */
   public int getCompressionRatio() {
     return compressionMultiplier*minCompressionRatio;
   }
   
   /** Retrieves compressed value associated with given index.
    *  <p>
    *  
    *  Indexes range from <tt>0...nElements-1</tt>, with order corresponding
    *  to the order of the underlying, uncompressed, data in the
    *  original sequence (<tt>0</tt> for the oldest,
    *  <tt>nElements-1</tt> for the newest).
    *
    * @param iCompressed index of compressed value to be retrieved
    * @return value of the compressed sequence at the given index, or
    * <tt>Double.NaN</tt> if no data has been placed into the
    * compressed buffer at the given index (occurs only during startup).
    *
    */ 
   public double getCompressedValue(int iCompressed) {
     if (iCompressed < 0 || iCompressed >= nElements)
       throw new IllegalArgumentException("iCompressed=" + iCompressed +
          "; iCompressed must be >= 0 and < nElements (which is " + nElements + ")");
     double result = Double.NaN;
     int index = getFirstSampleIndexUnchecked(iCompressed);
     if (index != NOSUCHINDEX) {
       index %= nElements;
       result = buffer[index];
     }
     return result;
   }

   /**
    * Retrieves the "in progress", or incomplete compressed value. 
    * <p>
    * 
    * This value reflects any recent updates not yet incorporated into
    * a compressed value. For example, if the buffer
    * is currently compressing at a 4 to 1 ratio, and the total number of updates
    * were not a multiple of 4, then a compressed value
    * for either the most recently seen, 1, 2, or 3 points would be returned.
    * 
    * <p>
    *
    * Note that if no incomplete compressed value is available,
    * <tt>Double.NaN</tt> will be returned.
    *
    * @return value of the incomplete or "in progress" compressed element.
    * 
    *
    */
   public double getIncompleteCompressedValue() {
     double result = compressionMethod.getCompressedValue();
     return result;
   }

   /**
    * Returns the compressing circular buffer to it's initial state,
    * discarding all previously seen data.
    * 
    */
   public void reset() {
     iElement = -1;  
     compressionMultiplier = 1;
     compressionMethod.reset();
   }
}
