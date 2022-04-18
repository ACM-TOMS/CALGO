package com.googlecode.compressingcircularbuffers;
/** Compresses by replacing adjacent values with a random sample
 ** thereof. */
public class RandomSamplingCompressionMethod implements CompressionMethod {
  int count = 0;
  private double randomSample = Double.NaN;
  public void update(double value) {
    count++;
        // All values seen after reset have equal chance of being
        // selected. For example, after 4 updates (count=4) we have,
        // from first to last, equal selection probabilities of
        // 1*(1/2)*(2/3)*(3/4), (1/2)*(2/3)*(3/4), (1/3)*(3/4), and 1/4.
    if (count == 1 || Math.random() < 1./count)
      randomSample = value;
  }
  public void reset() {
    count = 0;
    randomSample = Double.NaN;
  }
  public double getCompressedValue() {
    return randomSample;
  }
  public double mergeValues(double older, double newer) {
    double result = (Math.random() < 0.5) ? older : newer;
    return result;
  }
  public int getCount() { return count; }
}
