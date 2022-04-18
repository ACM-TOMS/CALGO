package com.googlecode.compressingcircularbuffers;

/** Compresses by aggregating adjacent samples into averages. */
public class AveragingCompressionMethod implements CompressionMethod {
  private double sum = 0;
  private int count = 0;
  public void update(double value) {
    sum += value;
    count++;
  }
  public void reset() {
    sum = 0;
    count = 0;
  }
  public double getCompressedValue() {
    double result = Double.NaN;
    if (count > 0)
      result = sum/count;
    return result;
  }
  public double mergeValues(double older, double newer) {
    double result = 0.5*(older + newer);
    return result;
  }
  public int getCount() { return count; }
}
