package com.googlecode.compressingcircularbuffers;
/**
 * A compression method defines how <tt>2<sup>k</sup></tt> adjacent samples get
 * transformed into a single compressed value.
 * <p>
 *
 * Rather than implement this interface, most applications can use one of
 * the predefined compression methods listed below.
 *
 * @see ThinningCompressionMethod ThinningCompressionMethod
 * @see AveragingCompressionMethod AveragingCompressionMethod
 * @see SummingCompressionMethod SummingCompressionMethod
 * @see MinimizingCompressionMethod MinimizingCompressionMethod
 * @see MaximizingCompressionMethod MaximizingCompressionMethod
 * @see RandomSamplingCompressionMethod RandomSamplingCompressionMethod
 * 
 *
 */
public interface CompressionMethod {
/**
 *  Modifies compressed value to reflect one additional sample value.
 * 
 *  @param value sampled value used to update the compressed value.
 * 
 */
  void update(double value);
/** Combines two adjacent n-compressed values to form a single
 *  2n-compressed value.
 * 
 *  @param older the older of the two adjacent values.
 *  @param newer the newer of the two adjacent values.
 * 
 *  @return single value that combines the two values. For example, an
 *  averaging compression method would return <tt>(older+newer)/2</tt>, whereas a
 *  maximizing compression method would return
 *  <tt>max(older,newer)</tt>.
 * 
 */
  double mergeValues(double older, double newer);
  
/** Eliminates the impact of previously seen data from the compression
 *  method. Each time a new compressed value is formed, it is recorded within
 *  <tt>CompressionCircularBuffer</tt>'s buffer, and this method is
 *  then called to make way for a fresh chunk of data.
 *
 */
  void reset();
/**
 *  Compressed value representing all updates since last reset of this
 *  compression method.
 * 
 *  @return current compressed value associated with this compression
 *  method.
 *
 */
  double getCompressedValue();
/**
 *  Returns number of samples (updates) in current compressed value.
 * 
 *  @return number of samples since this compression method was last
 *  reset.
 * 
 */
  int getCount();
}
