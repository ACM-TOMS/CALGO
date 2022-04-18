package com.googlecode.compressingcircularbuffers;
/** Compresses by replacing adjacent values with minimum thereof. */
public class MinimizingCompressionMethod implements CompressionMethod {
  int count = 0;
  double min = Double.NaN;
  public void update(double value) {
    if (0 == count || value < min)
      min = value;          
    count++;
  }
  public void reset() {
    min = Double.NaN;
    count = 0;
  }
  public double getCompressedValue() {
    return min;
  }
  public double mergeValues(double older, double newer) {
    double result = Math.min(older, newer);
    return result;
  }
  public int getCount() { return count; }
}
