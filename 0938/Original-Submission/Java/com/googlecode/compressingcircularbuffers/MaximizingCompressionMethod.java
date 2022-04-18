package com.googlecode.compressingcircularbuffers;
/** Compresses by replacing adjacent values with maximum thereof. */
public class MaximizingCompressionMethod implements CompressionMethod {
  int count = 0;
  double max = Double.NaN;
  public void update(double value) {
    if (0 == count || value > max)
      max = value;
    count++;
  }
  public void reset() {
    max = Double.NaN;
    count = 0;
  }
  public double getCompressedValue() {
    return max;
  }
  public double mergeValues(double older, double newer) {
    double result = Math.max(older, newer);
    return result;
  }
  public int getCount() { return count; }
}
